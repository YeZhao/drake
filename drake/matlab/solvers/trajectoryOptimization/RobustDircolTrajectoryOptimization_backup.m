classdef RobustDircolTrajectoryOptimization_backup < DirectTrajectoryOptimization
  % Direct colocation approach
  % Over each interval, f(x(k),u(k)) and f(x(k+1),u(k+1)) are evaluated to
  % determine d/dt x(k) and d/dt x(k+1). A cubic spline is fit over the
  % interval x and d/dt x at the end points.
  % x(k+.5) and d/dt x(k+.5) are determined based on this spline.
  % Then, the dynamics constraint is:
  % d/dt x(k+.5) = f(x(k+.5),.5*u(k) + .5*u(k+1))
  %
  %  integrated cost is: .5*h(1)*g(x(1),u(1)) + .5*h(N-1)*g(x(N),u(N)) +
  %                   sum((.5*h(i)+.5*h(i-1))*g(x(i),u(i))
  %  more simply stated, integrated as a zoh with half of each time
  %  interval on either side of the knot point
  % this might be the wrong thing for the cost function...
  properties
      nX
    nU
    nW
    
    D % Disturbance ellipsoid matrix w'*D^-1*w <= 1
    E0 % Initial disturbed state at t=0
    
    Q % LQR state cost matrix
    R % LQR input cost matrix
    Qf% LQR terminal cost matrix
    
    Qr %Robust cost matrix
    Rr %Robust cost matrix
    Qrf %Robust cost matrix
    
    xtraj_nom % State trajectory of nominal model (non-perturbed) 
    utraj_nom % Control trajectory of nominal model (non-perturbed)
    K_nom % LQR gain matrix of nominal model (non-perturbed)
    
    xtraj_perturbed_array % State trajectory array of perturbed model
    utraj_perturbed_array % Control trajectory array of perturbed model
    
    Pc %Projection onto constrained subspace of state vector
    
    %Stuff to cache so we don't have to recompute LQR controller
    z_handle
    K_handle
    A_handle
    B_handle
    G_handle
    dK_handle
    dA_handle
    dB_handle
    dG_handle
    E_handle
    dE_handle
    
    %Stuff for robustifying state constraints
    constr_xinds
    dUdE
    dxcdv
    dxc
  end
  
  methods
    function obj = RobustDircolTrajectoryOptimization_backup(plant,N,duration,varargin)
      obj = obj@DirectTrajectoryOptimization(plant,N,duration,varargin{:});
      
      obj.nX = plant.getNumStates;
      obj.nU = plant.getNumInputs;
      obj.nW = plant.getNumDisturbances;
      obj.D = D;
      obj.E0 = E0;
      obj.Q = Q;
      obj.R = R;
      obj.Qf = Qf;
      
      obj.z_handle = SharedDataHandle(0);
      obj.K_handle = SharedDataHandle(0);
      obj.A_handle = SharedDataHandle(0);
      obj.B_handle = SharedDataHandle(0);
      obj.G_handle = SharedDataHandle(0);
      obj.dK_handle = SharedDataHandle(0);
      obj.dA_handle = SharedDataHandle(0);
      obj.dB_handle = SharedDataHandle(0);
      obj.dG_handle = SharedDataHandle(0);
      obj.E_handle = SharedDataHandle(0);
      obj.dE_handle = SharedDataHandle(0);
    end
    
    function obj = setupVariables(obj, N)
      % Assumes that there are N-1 time steps
      % N corresponding state variables
      % and N-1 corresponding input variables
      %
      % Generates num_vars total number of decision variables
      %   h_inds (N-1) x 1 indices for timesteps h so that z(h_inds(i)) = h(i)
      %   x_inds N x n indices for state
      %   u_inds (N-1) x m indices for input
      %
      % @param N number of knot points
      nH = N-1;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();

      num_vars = nH + N*nX + (N-1)*nU;
      obj.h_inds = (1:nH)';
      obj.x_inds = reshape(nH + (1:nX*N),nX,N);
      obj.u_inds = reshape(nH + nX*N + (1:nU*(N-1)),nU,N-1);

      obj.N = N;
      x_names = cell(num_vars,1);
      for i = 1:(N-1)
        x_names{i} = sprintf('h[%d]',i);
        for j = 1:nX
          x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
        end
        for j = 1:nU
          x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
        end
      end
      for j = 1:nX
          x_names{nH+(N-1)*nX+j}=sprintf('x%d[%d]',j,N);
      end

      obj = obj.addDecisionVariable(num_vars,x_names);
      
      % Ensure that all h values are non-negative
      obj = obj.addConstraint(BoundingBoxConstraint(zeros(N-1,1),inf(N-1,1)),obj.h_inds);
      
      % create constraints for dynamics and add them
      obj = obj.addDynamicConstraints();
      
      % add control inputs as bounding box constraints
      if any(~isinf(obj.plant.umin)) || any(~isinf(obj.plant.umax))
          control_limit = BoundingBoxConstraint(repmat(obj.plant.umin,N-1,1),repmat(obj.plant.umax,N-1,1));
          obj = obj.addConstraint(control_limit,obj.u_inds(:));
      end
      
      obj.z_handle.data = zeros((N-1)*(1+nX+nU)+nX,1);%randn((N-1)*(1+nX+nU)+nX,1);
      
    end
    
    function z0 = getInitialVars(obj,t_init,traj_init)
        % evaluates the initial trajectories at the sampled times and
        % constructs the nominal z0.
        if isscalar(t_init)
            t_init = linspace(0,t_init,obj.N);
        elseif length(t_init) ~= obj.N
            error('The initial sample times must have the same length as property N')
        end
        z0 = zeros(obj.num_vars,1);
        z0(obj.h_inds) = diff(t_init);
        
        if nargin<3, traj_init = struct(); end
        
        nU = getNumInputs(obj.plant);
        if isfield(traj_init,'u')
            z0(obj.u_inds) = traj_init.u.eval(t_init(1:end-1));
        else
            z0(obj.u_inds) = zeros(nU,obj.N-1);%0.01*randn(nU,obj.N-1);
        end
        
        if isfield(traj_init,'x')
            z0(obj.x_inds) = traj_init.x.eval(t_init);
        else
            if nU>0
                if ~isfield(traj_init,'u')
                    traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nU,obj.N-1))),getInputFrame(obj.plant));
                end
                
                % todo: if x0 and xf are equality constrained, then initialize with
                % a straight line from x0 to xf (this was the previous behavior)
                
                %simulate
                sys_ol = cascade(traj_init.u,obj.plant);
            else
                sys_ol = obj.plant;
            end
            
            if ~isfield(traj_init,'x0')
                [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)]);
            else
                [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)],traj_init.x0);
            end
            
            z0(obj.x_inds) = x_sim.eval(t_init);
        end
    end
    
    function obj = addDynamicConstraints(obj)
      N = obj.N;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      constraints = cell(N-1,1);
      dyn_inds = cell(N-1,1);
      
      
      n_vars = 2*nX + 2*nU + 1;
      cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.constraint_fun);
      cnstr = setName(cnstr,'collocation');

      % create shared data functions to calculate dynamics at the knot
      % points
      shared_data_index = obj.getNumSharedDataFunctions;
      for i=1:obj.N,
        obj = obj.addSharedDataFunction(@obj.dynamics_data,{obj.x_inds(:,i);obj.u_inds(:,i)});
      end
      
      for i=1:obj.N-1,
        dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
        constraints{i} = cnstr;
        
        obj = obj.addConstraint(constraints{i}, dyn_inds{i},[shared_data_index+i;shared_data_index+i+1]);
      end
    end
    
    function [f,df] = constraint_fun(obj,h,x0,x1,u0,u1,data0,data1)
      % calculate xdot at knot points
      %  [xdot0,dxdot0] = obj.plant.dynamics(0,x0,u0);
      %  [xdot1,dxdot1] = obj.plant.dynamics(0,x1,u1);
      
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      
      % use the shared data objects for the dynamics at the knot points
      xdot0 = data0.xdot;
      dxdot0 = data0.dxdot;
      xdot1 = data1.xdot;
      dxdot1 = data1.dxdot;
      
      % cubic interpolation to get xcol and xdotcol, as well as
      % derivatives
      xcol = .5*(x0+x1) + h/8*(xdot0-xdot1);
      dxcol = [1/8*(xdot0-xdot1) (.5*eye(nX) + h/8*dxdot0(:,2:1+nX)) ...
        (.5*eye(nX) - h/8*dxdot1(:,2:1+nX)) h/8*dxdot0(:,nX+2:1+nX+nU) -h/8*dxdot1(:,nX+2:1+nX+nU)];
      xdotcol = -1.5*(x0-x1)/h - .25*(xdot0+xdot1);
      dxdotcol = [1.5*(x0-x1)/h^2 (-1.5*eye(nX)/h - .25*dxdot0(:,2:1+nX)) ...
        (1.5*eye(nX)/h - .25*dxdot1(:,2:1+nX)) -.25*dxdot0(:,nX+2:1+nX+nU) -.25*dxdot1(:,nX+2:1+nX+nU)];
      
      % evaluate xdot at xcol, using foh on control input
      [g,dgdxcol] = obj.plant.dynamics(0,xcol,.5*(u0+u1));
      dg = dgdxcol(:,2:1+nX)*dxcol + [zeros(nX,1+2*nX) .5*dgdxcol(:,2+nX:1+nX+nU) .5*dgdxcol(:,2+nX:1+nX+nU)];
      
      % constrait is the difference between the two
      f = xdotcol - g;
      df = dxdotcol - dg;
    end
    
    
    function data = dynamics_data(obj,x,u)
      [data.xdot,data.dxdot] = obj.plant.dynamics(0,x,u);
    end
    
    function obj = addRunningCost(obj,running_cost_function)
      % Adds an integrated cost to all time steps, which is
      % numerical implementation specific (thus abstract)
      % this cost is assumed to be time-invariant
      % @param running_cost_function a function handle
      %  of the form running_cost_function(dt,x,u)
      % This implementation assumes a ZOH, but where the values of
      % x(i),u(i) are held over an interval spanned by .5(dt(i-1) + dt(i))
      
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
            
      running_cost_end = FunctionHandleObjective(1+nX+nU,@(h,x,u) obj.running_fun_end(running_cost_function,h,x,u));
      running_cost_mid = FunctionHandleObjective(2+nX+nU,@(h0,h1,x,u) obj.running_fun_mid(running_cost_function,h0,h1,x,u));

      obj = obj.addCost(running_cost_end,{obj.h_inds(1);obj.x_inds(:,1);obj.u_inds(:,1)});
      
      for i=2:obj.N-1,
        obj = obj.addCost(running_cost_mid,{obj.h_inds(i-1);obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)});
      end
      
      obj = obj.addCost(running_cost_end,{obj.h_inds(end);obj.x_inds(:,end);obj.u_inds(:,end)});
    end
    
    function [utraj,xtraj] = reconstructInputTrajectory(obj,z)
      % Interpolate between knot points to reconstruct a trajectory using
      % the hermite spline
      t = [0; cumsum(z(obj.h_inds))];
      u = reshape(z(obj.u_inds),[],obj.N);
      utraj = PPTrajectory(foh(t,u));
      utraj = utraj.setOutputFrame(obj.plant.getInputFrame);
    end
    
    function xtraj = reconstructStateTrajectory(obj,z)
      % Interpolate between knot points to reconstruct a trajectory using
      % the hermite spline
      t = [0; cumsum(z(obj.h_inds))];
      u = reshape(z(obj.u_inds),[],obj.N);

      x = reshape(z(obj.x_inds),[],obj.N);
      xdot = zeros(size(x,1),obj.N);
      for i=1:obj.N,
        xdot(:,i) = obj.plant.dynamics(t(i),x(:,i),u(:,i));
      end
      xtraj = PPTrajectory(pchipDeriv(t,x,xdot));
      xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
    end
  end
  
  methods (Access = protected)
    function [f,df] = running_fun_end(obj,running_handle,h,x,u)
      [f,dg] = running_handle(.5*h,x,u);
      
      df = [.5*dg(:,1) dg(:,2:end)];
    end
    
    function [f,df] = running_fun_mid(obj,running_handle,h0,h1,x,u)
      [f,dg] = running_handle(.5*(h0+h1),x,u);
      
      df = [.5*dg(:,1) .5*dg(:,1) dg(:,2:end)];
    end
    
  end
end