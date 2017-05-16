classdef RobustDircolTrajectoryOptimization < DirectTrajectoryOptimization
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
      
      Qr %Robust cost matrix
      Rr %Robust cost matrix
      Qrf %Robust cost matrix
      
      xtraj_nom % State trajectory of nominal model (non-perturbed)
      utraj_nom % Control trajectory of nominal model (non-perturbed)
      K_nom % LQR gain matrix of nominal model (non-perturbed)
      
      xtraj_perturbed_array % State trajectory array of perturbed model
      utraj_perturbed_array % Control trajectory array of perturbed model
      
  end
  
  methods
    function obj = RobustDircolTrajectoryOptimization(plant,N,duration,varargin)
      obj = obj@DirectTrajectoryOptimization(plant,N,duration,varargin{:});
      
      obj.nX = plant.getNumStates;
      obj.nU = plant.getNumInputs;
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
    
    function obj = addRobustRunningCost(obj,utraj,xtraj,K,Qr,Qrf,Rr)
        nX = obj.nX;
        nU = obj.nU;
        N = obj.N;
        
        obj.Qr = Qr;
        obj.Rr = Rr;
        obj.Qrf = Qrf;
        obj.xtraj_nom = xtraj;
        obj.utraj_nom = utraj;
        obj.K_nom = K;
        
        dim = N-1 + N*nX + (N-1)*nU;
        cost = FunctionHandleObjective(dim,@obj.robust_cost_sampled,1);
        cost.grad_method = 'user';
        obj = obj.addCost(cost, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds(:,1:end-1)],[],1); obj.x_inds(:,end)});
    end
    
    function obj = addRobustAverageRunningCost(obj,utrajArray,xtrajArray,Qr,Qrf,Rr)
        
        obj.Qr = Qr;
        obj.Rr = Rr;
        obj.Qrf = Qrf;
        obj.xtraj_perturbed_array = xtrajArray;
        obj.utraj_perturbed_array = utrajArray;
        
        dim = obj.N-1 + obj.N*obj.nX + (obj.N-1)*obj.nU;
        cost = FunctionHandleObjective(dim,@obj.robust_cost_average,1);
        cost.grad_method = 'user';
        obj = obj.addCost(cost, {reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds(:,1:end-1)],[],1); obj.x_inds(:,end)});
    end
    
    function [c, dc] = robust_cost_sampled(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        N = obj.N;
                
        %decompose the decision variable vector y        
        for k = 1:N-1
            xtraj(:,k) = y((k-1)*(1+nX+nU)+1+(1:nX));
            utraj(:,k) = y((k-1)*(1+nX+nU)+1+nX+(1:nU));
        end
        xtraj(:,N) = xf;
        
        c = 0;
        dc = zeros(1,(N-1)*(1+nX+nU)+nX);
        for k = 1:(N-1)
            % method 1: do not enforce the LQR feedback
            %c = c + .5*(xtraj(:,k) - obj.xtraj_nom(:,k))'*obj.Qr*(xtraj(:,k) - obj.xtraj_nom(:,k)) + ...
            %    .5*(utraj(:,k) - obj.utraj_nom(:,k))'*obj.Rr*(utraj(:,k) - obj.utraj_nom(:,k));
            
            % method 2: enforce the LQR feedback
            c = c + .5*(xtraj(:,k) - obj.xtraj_nom(:,k))'*(obj.Qr + obj.K_nom(:,:,k)'*obj.Rr*obj.K_nom(:,:,k)) ...
                *(xtraj(:,k) - obj.xtraj_nom(:,k));
            
            dc((k-1)*(1+nX+nU)+1) = 0;
            % method 1: do not enforce the LQR feedback
            %dc((k-1)*(1+nX+nU)+1+(1:nX)) = (xtraj(:,k) - obj.xtraj_nom(:,k))'*obj.Qr;
            %dc((k-1)*(1+nX+nU)+1+nX+(1:nU)) = (utraj(:,k) - obj.utraj_nom(:,k))'*obj.Rr;
            
            % method 2: enforce the LQR feedback
            dc((k-1)*(1+nX+nU)+1+(1:nX)) = (xtraj(:,k) - obj.xtraj_nom(:,k))'*(obj.Qr + obj.K_nom(:,:,k)'*obj.Rr*obj.K_nom(:,:,k));
            dc((k-1)*(1+nX+nU)+1+nX+(1:nU)) = zeros(1,nU);
        end
        c = c + .5*(xtraj(:,N) - obj.xtraj_nom(:,N))'*obj.Qrf*(xtraj(:,N) - obj.xtraj_nom(:,N));
        dc((N-1)*(1+nX+nU)+(1:nX)) = (xtraj(:,N) - obj.xtraj_nom(:,N))'*obj.Qrf;
        
%         persistent count_robust;
%         if isempty(count_robust)
%             count_robust = 1;
%         else
%             count_robust = count_robust + 1;
%         end
%         count_robust
%         c
    end
    
    function [c, dc] = robust_cost_average(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        N = obj.N;
                
        %decompose the decision variable vector y        
        for k = 1:N-1
            xtraj(:,k) = y((k-1)*(1+nX+nU)+1+(1:nX));
            utraj(:,k) = y((k-1)*(1+nX+nU)+1+nX+(1:nU));
        end
        xtraj(:,N) = xf;
        
        c = 0;
        dc = zeros(1,(N-1)*(1+nX+nU)+nX);
        Numoftraj = length(obj.xtraj_perturbed_array(1,1,:));
        for j = 1:Numoftraj
            for k = 1:(N-1)
                c = c + .5*(xtraj(:,k) - obj.xtraj_perturbed_array(:,k,j))'*obj.Qr*(xtraj(:,k) - obj.xtraj_perturbed_array(:,k,j)) + ...
                   .5*(utraj(:,k) - obj.utraj_perturbed_array(:,k,j))'*obj.Rr*(utraj(:,k) - obj.utraj_perturbed_array(:,k,j));

                dc((k-1)*(1+nX+nU)+1) = 0;
                dc((k-1)*(1+nX+nU)+1+(1:nX)) = dc((k-1)*(1+nX+nU)+1+(1:nX)) + (xtraj(:,k) - obj.xtraj_perturbed_array(:,k,j))'*obj.Qr;
                dc((k-1)*(1+nX+nU)+1+nX+(1:nU)) = dc((k-1)*(1+nX+nU)+1+nX+(1:nU)) + (utraj(:,k) - obj.utraj_perturbed_array(:,k,j))'*obj.Rr;                
            end
            c = c + .5*(xtraj(:,N) - obj.xtraj_perturbed_array(:,N,j))'*obj.Qrf*(xtraj(:,N) - obj.xtraj_perturbed_array(:,N,j));
            dc((N-1)*(1+nX+nU)+(1:nX)) = dc((N-1)*(1+nX+nU)+(1:nX)) + (xtraj(:,N) - obj.xtraj_perturbed_array(:,N,j))'*obj.Qrf;
        end
        c = c/Numoftraj;
        dc = dc/Numoftraj;
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