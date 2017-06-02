classdef IDQP < DrakeSystem
  methods
  function obj = IDQP(r,xtraj,utraj,ctraj,options)
    % @param r rigid body manipulator instance
    % @param options structure for specifying objective weights, slack
    % bounds, etc.
    
    if nargin>4
      typecheck(options,'struct');
    else
      options = struct();
    end

    input_frame = r.getStateFrame();
    output_frame = r.getInputFrame();
    
    obj = obj@DrakeSystem(0,0,r.getNumStates,r.getNumInputs,true,true);

    obj = setInputFrame(obj,input_frame);
    obj = setOutputFrame(obj,output_frame);

    obj.robot = r;
    obj.numq = getNumPositions(r);
    obj.numv = getNumVelocities(r);
    
    if isfield(options,'dt')
      % controller update rate
      typecheck(options.dt,'double');
      sizecheck(options.dt,[1 1]);
      dt = options.dt;
    else
      dt = 0.001;
    end
    obj = setSampleTime(obj,[dt;0]); % sets controller update rate

    if isfield(options,'use_bullet')
      obj.use_bullet = options.use_bullet;
    else
      obj.use_bullet = false;
    end

    % weight for the desired qddot objective term
    if isfield(options,'w_qdd')
      typecheck(options.w_qdd,'double');
      sizecheck(options.w_qdd,[obj.numq 1]); % assume diagonal cost
      obj.w_qdd = options.w_qdd;
    else
      obj.w_qdd = 1*ones(obj.numq,1);
    end

    % weight for slack vars
    if isfield(options,'w_slack')
      typecheck(options.w_slack,'double');
      sizecheck(options.w_slack,1);
      obj.w_slack = options.w_slack;
    else
      obj.w_slack = 0.001;
    end

    % gain for support acceleration constraint: accel=-Kp_accel*vel
    if isfield(options,'Kp_accel')
      typecheck(options.Kp_accel,'double');
      sizecheck(options.Kp_accel,1);
      obj.Kp_accel = options.Kp_accel;
    else
      obj.Kp_accel = 0.0; % default desired acceleration=0
    end

    % hard bound on slack variable values
    if isfield(options,'slack_limit')
      typecheck(options.slack_limit,'double');
      sizecheck(options.slack_limit,1);
      obj.slack_limit = options.slack_limit;
    else
      obj.slack_limit = 0.1;
    end

    if isfield(options,'solver')
      % 0: fastqp, fallback to gurobi barrier (default)
      % 1: gurobi primal simplex with active sets
      typecheck(options.solver,'double');
      sizecheck(options.solver,1);
      assert(options.solver==0 || options.solver==1);
    else
      options.solver = 0;
    end
    obj.solver = options.solver;

    obj.gurobi_options.outputflag = 0; % not verbose
    if options.solver==0
      obj.gurobi_options.method = 2; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    else
      obj.gurobi_options.method = 0; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
    end
    obj.gurobi_options.presolve = 0;
    % obj.gurobi_options.prepasses = 1;

    if obj.gurobi_options.method == 2
      obj.gurobi_options.bariterlimit = 20; % iteration limit
      obj.gurobi_options.barhomogeneous = 0; % 0 off, 1 on
      obj.gurobi_options.barconvtol = 5e-4;
    end

    obj.using_flat_terrain = true;
    [obj.jlmin, obj.jlmax] = getJointLimits(r);
    
    obj.xtraj = xtraj;
    obj.utraj = utraj;
    obj.ctraj = ctraj;
  end

  function y=output(obj,t,~,x)

    r = obj.robot;
    nq = obj.numq;
    nu = r.getNumInputs;
    q = x(1:nq);
    qd = x(nq+(1:nq));
    
    kinsol = doKinematics(r,q,false,true,qd);

    dim = 2; % 3D
    nd = 2; % for friction cone approx, hard coded for now

    [H,C,B] = manipulatorDynamics(r,q,qd);

%     [xcom,Jcom] = getCOM(r,kinsol);
%     Jcomdot_times_v = centerOfMassJacobianDotTimesV(r, kinsol);
    
    [phi,~,JB] = contactConstraintsBV(r,kinsol,false,struct('terrain_only',~obj.use_bullet));

    active = phi<=0.001;
    nc = sum(active);
    Dbar = vertcat(JB{active})';
    
%     terrain_pts = getTerrainContactPoints(r);
%     [~,Jp] = terrainContactPositions(r,kinsol,terrain_pts,true);
%     Jpdot_times_v = terrainContactJacobianDotTimesV(r, kinsol, terrain_pts);
%     Jp = sparse(Jp);
   
    neps = nc*dim;

    %----------------------------------------------------------------------
    % Build handy index matrices ------------------------------------------

    nf = nc*nd; % number of contact force variables
    nparams = nq+nu+nf;
    Iqdd = zeros(nq,nparams); Iqdd(:,1:nq) = eye(nq);
    Iu = zeros(nu,nparams); Iu(:,nq+(1:nu)) = eye(nu);
    Ibeta = zeros(nf,nparams); Ibeta(:,nq+nu+(1:nf)) = eye(nf);

    %----------------------------------------------------------------------
    % Set up problem constraints ------------------------------------------

    lb = [-inf(1,nq) r.umin'  zeros(1,nf)]'; % qddot/u/contact+constraint forces/slack vars
    ub = [inf(1,nq) r.umax' inf(1,nf)]';

    Aeq_ = cell(1,4);
    beq_ = cell(1,4);

    % constrained dynamics
    if nc>0
      Aeq_{1} = H*Iqdd - B*Iu - Dbar*Ibeta;
    else
      Aeq_{1} = H*Iqdd - B*Iu;
    end
    beq_{1} = -C;

    % linear equality constraints: Aeq*alpha = beq
    Aeq = sparse(vertcat(Aeq_{:}));
    beq = vertcat(beq_{:});

%     % linear inequality constraints: Ain*alpha <= bin
%     Ain = sparse(vertcat(Ain_{:}));
%     bin = vertcat(bin_{:});
%     Ain = Ain(bin~=inf,:);
%     bin = bin(bin~=inf);
    Ain = [];
    bin = [];

    Kp_qdd = 500*eye(nq);
%     Kd_qdd = 1*eye(nq);
    Kd_qdd = 1*sqrt(Kp_qdd);

    idx = [r.findJointId('left_ankle');
          r.findJointId('right_ankle')]-1;
    Kp_qdd(idx,idx) = 0;
    Kd_qdd(idx,idx) = 0;

    w_q = obj.w_qdd;
    w_q(idx)=0;
    w_q(1)=5*w_q(1);
    
    x_des = obj.xtraj.eval(t);
    qddot_des = Kp_qdd*(x_des(1:nq) - q) + Kd_qdd*(x_des(nq+(1:nq))-qd);

    u_des = obj.utraj.eval(t);
    w_u = .01;
    
    c_des = obj.ctraj.eval(t);
    c_des = c_des(active);
    beta_des = kron(eye(nc),[1;1]) * c_des;
    w_beta = 1e-4;

    %----------------------------------------------------------------------
    % QP cost function ----------------------------------------------------
    %
    Hqp = zeros(nparams);
    Hqp(1:nq,1:nq) = diag(w_q);
    Hqp(nq+(1:nu),nq+(1:nu)) = w_u*eye(nu);
    Hqp(nq+nu+(1:nf),nq+nu+(1:nf)) =  w_beta*eye(nf);
    fqp = -(w_q.*qddot_des)'*Iqdd;
    fqp = fqp - (w_u.*u_des)'*Iu;
    fqp = fqp - (w_beta.*beta_des)'*Ibeta;

    %----------------------------------------------------------------------
    % Solve QP ------------------------------------------------------------

    REG = 1e-8;

    IR = eye(nparams);
    lbind = lb>-999;  ubind = ub<999;  % 1e3 was used like inf above... right?
    Ain_fqp = full([Ain; -IR(lbind,:); IR(ubind,:)]);
    bin_fqp = [bin; -lb(lbind); ub(ubind)];

%     % call fastQPmex first
%     QblkDiag = {Hqp(1:nq,1:nq) + REG*eye(nq), ...
%                 REG*ones(nu,1), ...
%                 REG*ones(nf,1), ...
%                 obj.w_slack*ones(neps,1) + REG*ones(neps,1)};
%     Aeq_fqp = full(Aeq);
%     % NOTE: model.obj is 2* f for fastQP!!!
%     [alpha,info_fqp] = fastQPmex(QblkDiag,fqp,Ain_fqp,bin_fqp,Aeq_fqp,beq);
% 
    if 1% info_fqp<0
      % then call gurobi
      % disp('QPController: failed over to gurobi');
      model.Q = sparse(Hqp + REG*eye(nparams));
      model.A = [Aeq; Ain];
      model.rhs = [beq; bin];
      model.sense = [obj.eq_array(1:length(beq)); obj.ineq_array(1:length(bin))];
      model.lb = lb;
      model.ub = ub;

      model.obj = fqp;
      if obj.gurobi_options.method==2
        % see drake/algorithms/QuadraticProgram.m solveWGUROBI
        model.Q = .5*model.Q;
      end

      if (any(any(isnan(model.Q))) || any(isnan(model.obj)) || any(any(isnan(model.A))) || any(isnan(model.rhs)) || any(isnan(model.lb)) || any(isnan(model.ub)))
        keyboard;
      end

%         qp_tic = tic;
      result = gurobi(model,obj.gurobi_options);
%         qp_toc = toc(qp_tic);
%         fprintf('QP solve: %2.4f\n',qp_toc);

      alpha = result.x;
    end
    
    
    
    y = alpha(nq+(1:nu));
    
      x_des - x
    u_des-y
    

  end

  end
  
  properties (SetAccess=private)
    robot; % to be controlled
    numq;
    numv;
    controller_data; % shared data handle that holds S, h, foot trajectories, etc.
    W_kdot; % angular momentum cost term weight matrix
    w_qdd; % qdd objective function weight vector
    w_grf; % scalar ground reaction force weight
    w_slack; % scalar slack var weight
    slack_limit; % maximum absolute magnitude of acceleration slack variable values
    Kp_ang; % proportunal gain for angular momentum feedback
    Kp_accel; % gain for support acceleration constraint: accel=-Kp_accel*vel
    gurobi_options = struct();
    solver=0;
    lc;
    eq_array = repmat('=',100,1); % so we can avoid using repmat in the loop
    ineq_array = repmat('<',100,1); % so we can avoid using repmat in the loop
    use_bullet;
    using_flat_terrain; % true if using DRCFlatTerrain
    jlmin;
    jlmax;
    output_qdd = false;
    body_accel_input_weights; % array of doubles, negative values signal constraints
    n_body_accel_inputs;
    body_accel_bounds;
    n_body_accel_bounds;
    xtraj;
    utraj;
    ctraj;
    btraj;
  end
end
