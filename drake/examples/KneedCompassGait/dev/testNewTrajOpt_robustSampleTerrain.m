function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testNewTrajOpt_robustSampleTerrain(xtraj,utraj,ltraj,ljltraj)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');

paramstd = 1/5; % Standard deviation of the parameter value percent error
SampleNum = 1; % number of sampled terrain height
% perturb model parameters
paramerr = [];
for i = 1:SampleNum
    paramerr(i) = 0;%randn(1,1)*paramstd;
    if (paramerr(i) > 0.15)
        paramerr(i) = 0.15;
    elseif (paramerr(i) < -0.15)
        paramerr(i) = -0.15;
    end
    options(i).terrain = RigidBodyStepTerrainVaryingHeight(paramerr(i));
    options(i).floating = true;
    options(i).ignore_self_collisions = true;
    
    %p_perturb(i) = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options(i));
    p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options(i));
    
    % trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);
end

v = p.constructVisualizer;

%todo: add joint limits, periodicity constraint

N = 10;
T = 5;
T0 = 5;

% periodic constraint
R_periodic = zeros(p.getNumStates,2*p.getNumStates);
R_periodic(2,2) = 1; %z
R_periodic(3,3) = 1; %pitch-hip w/symmetry
R_periodic(3,5) = 1; %pitch-hip w/symmetry
R_periodic(4,6) = 1; %knee w/symmetry
R_periodic(6,4) = 1; %knee w/symmetry
R_periodic(5,5) = -1; %hip w/symmetry

R_periodic(7,7) = 1; %x-vel
R_periodic(8,8) = 1; %z-vel
R_periodic(9,9) = 1; %pitch-hip w/symmetry
R_periodic(9,11) = 1; %pitch-hip w/symmetry
R_periodic(10,12) = 1; %knee w/symmetry
R_periodic(12,10) = 1; %knee w/symmetry
R_periodic(11,11) = -1; %hip w/symmetry

R_periodic(2:end,p.getNumStates+2:end) = -eye(p.getNumStates-1);

periodic_constraint = LinearConstraint(zeros(p.getNumStates,1),zeros(p.getNumStates,1),R_periodic);

% x0 = [0;0;1;zeros(15,1)];
% xf = [0;0;1;zeros(15,1)];
x0 = [0;1;zeros(10,1)];
xf = [.2;1;zeros(10,1)];

N2 = floor(N/2);

if nargin < 2
    %Try to come up with a reasonable trajectory
    %x1 = [.3;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
    x1 = [.3;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
    t_init = linspace(0,T0,N);
    %   traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
    traj_init.x = PPTrajectory(foh(t_init,[linspacevec(x0,x1,N2), linspacevec(x1,xf,N-N2)]));
    traj_init.u = PPTrajectory(foh(t_init,randn(3,N)));
    traj_init.l = PPTrajectory(foh(t_init,[repmat([1;zeros(7,1)],1,N2) repmat([zeros(4,1);1;zeros(3,1)],1,N-N2)]));
    traj_init.ljl = PPTrajectory(foh(t_init,zeros(p.getNumJointLimitConstraints,N)));
else
    t_init = xtraj.pp.breaks;
    traj_init.x = xtraj;
    traj_init.u = utraj;
    traj_init.l = ltraj;
    traj_init.ljl = ljltraj;
end
T_span = [1 T];


x0_min = [x0(1:5);-inf; -inf; 0; -inf(4,1)];
x0_max = [x0(1:5);inf;  inf; 0; inf(4,1)];
xf_min = [.4;-inf(11,1)];
xf_max = inf(12,1);

to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.compl_slack = .01;
to_options.lincompl_slack = .001;
to_options.jlcompl_slack = .01;
to_options.lambda_mult = p.getMass*9.81*T0/N;
to_options.lambda_jl_mult = T0/N;

traj_opt = ContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);
traj_opt = traj_opt.addStateConstraint(periodic_constraint,{[1 N]});

% % extra robust cost function
% traj_opt.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + traj_opt.nC,traj_opt.nC*(1+obj.nD),traj_opt.options.nlcc_mode,traj_opt.options.compl_slack);
% traj_opt.nonlincompl_slack_inds{i} = traj_opt.num_vars+1:traj_opt.num_vars + traj_opt.nonlincompl_constraints{i}.n_slack;
% traj_opt = traj_opt.addConstraint(traj_opt.nonlincompl_constraints{i},[traj_opt.x_inds(:,i+1);gamma_inds;lambda_inds]);


% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);
[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

v = p.constructVisualizer;
v.playback(xtraj,struct('slider',true));

disp('finish traj opt')

    function [f,df] = running_cost_fun(h,x,u)
        f = h*u'*u;
        df = [u'*u zeros(1,12) 2*h*u'];
    end

    function [f,df] = nonlincompl_fun(y)
        nq = p.getNumPositions;
        nv = p.getNumVelocities;
        x = y(1:nq+nv+traj_opt.nC);
        z = y(nq+nv+traj_opt.nC+1:end);
        gamma = x(nq+nv+1:end);
        
        q = x(1:nq);
        v = x(nq+1:nq+nv);
        
        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = p.contactConstraints(q,false,traj_opt.options.active_collision_options);
        
        f = zeros(traj_opt.nC*(1+traj_opt.nD),1);
        df = zeros(traj_opt.nC*(1+traj_opt.nD),nq+nv+traj_opt.nC*(2+traj_opt.nD));
        
        f(1:1+traj_opt.nD:end) = phi;
        df(1:1+traj_opt.nD:end,1:nq) = n;
        for j=1:traj_opt.nD,
            f(1+j:1+traj_opt.nD:end) = gamma+D{j}*v;
            df(1+j:1+traj_opt.nD:end,nq+nv+(1:traj_opt.nC)) = eye(size(D{j},1));  %d/dgamma
            df(1+j:1+traj_opt.nD:end,nq+(1:nv)) = D{j};%d/dv
            df(1+j:1+traj_opt.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
        end 
    end
end