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
      nW
      nq
      nv
      nu
      nx
      h
       
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
  end
  
  methods
    function obj = RobustDircolTrajectoryOptimization(plant,N,Q,R,Qf,duration,varargin)
      obj = obj@DirectTrajectoryOptimization(plant,N,duration,varargin{:});
      
      obj.nX = plant.getNumStates;
      obj.nU = plant.getNumInputs;
      obj.nW = plant.getNumDisturbances;
      
      obj.Q = Q;
      obj.R = R;
      obj.Qf = Qf;
    end
    
    function obj = addDynamicConstraints(obj)
      N = obj.N;
      nX = obj.plant.getNumStates();
      nU = obj.plant.getNumInputs();
      obj.nx = nX;
      obj.nq = nX/2;
      obj.nv = nX/2;
      obj.nu = nU;
      obj.nX = obj.nx;
      obj.nU = obj.nu;
      
      constraints = cell(N-1,1);
      dyn_inds = cell(N-1,1);
      
      n_vars = 2*nX + 2*nU + 1;
      cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.constraint_fun);
      cnstr = setName(cnstr,'collocation');

      % create shared data functions to calculate dynamics at the knot points
      shared_data_index = obj.getNumSharedDataFunctions;
      for i=1:obj.N,
       obj = obj.addSharedDataFunction(@obj.dynamics_data,{obj.x_inds(:,i);obj.u_inds(:,i)});
      end
      
      for i=1:obj.N-1,
         dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
         constraints{i} = cnstr;
         
         obj = obj.addConstraint(constraints{i}, dyn_inds{i},[shared_data_index+i;shared_data_index+i+1]);
      end
      
      obj = obj.addCost(FunctionHandleObjective(1,@(h_inds)getTimeStep(obj,h_inds),1),{obj.h_inds(1)});
      
      x_inds_stack = reshape(obj.x_inds,obj.N*nX,[]);
      u_inds_stack = reshape(obj.u_inds,obj.N*nU,[]);
      obj = obj.addCost(FunctionHandleObjective(obj.N*(nX+nU),@(x_inds,u_inds)robustVariancecost_scaled(obj,x_inds,u_inds),1),{x_inds_stack;u_inds_stack});
    end
    
    function [c,dc] = getTimeStep(obj, h)
        global timestep
        timestep = h;
        c = 0;
        dc = 0;
    end
    
    function [c,dc] = robustVariancecost_scaled(obj, x_full, u_full)
        global timestep
        
        x = reshape(x_full, obj.nx, obj.N);
        u = reshape(u_full, obj.nu, obj.N);
        nq = obj.plant.getNumPositions;
        nv = obj.plant.getNumVelocities;
        nx = nq+nv;
        nu = obj.plant.getNumInputs;
        h = timestep;
        
        %compute LQR gains
        % for k = 1:obj.N-1
        %     y((k-1)*(1+nx+nu)+1+(1:nx)) = x(:,k);
        %     y((k-1)*(1+nx+nu)+1+nx+(1:nu)) = u(:,k);
        % end
        % xf = x(:,obj.N);
        % [K] = deltaLQR(obj,y,xf);
                
        % sigma points
        Px = zeros(obj.nx,obj.nx,obj.N);
        Px(:,:,1) = obj.options.Px_coeff*eye(nx);
        Px_init = Px(:,:,1);
        w_noise = [];
        Pw = [];
        
        if strcmp(obj.plant.uncertainty_source,'physical_parameter_uncertainty')
            param_uncertainty = load('physical_param_pertubation.dat');
        elseif strcmp(obj.plant.uncertainty_source,'generate_new_parameter_uncertainty_set')
            paramstd = 1/5; % Standard deviation of the parameter value percent error
            for i = 1:2*(obj.nx)
                % Perturb original parameter estimates with random percentage error
                % normally distributed with standard dev = paramstd, and greater than -1
                param_uncertainty(i,:) = randn(1,10)*paramstd;
                while sum(param_uncertainty(i,:)<=-1)~=0
                    param_uncertainty(param_uncertainty(i,:)<-1) = randn(1,sum(param_uncertainty(i,:)<-1))*paramstd;
                end
            end
            save -ascii physical_param_pertubation.dat param_uncertainty
            w_noise = param_uncertainty;
        end
        
        % disturbance variance
        % currently only consider object horizontal 2D position and friction coefficient
        nw = 0;%size(Pw,1);
        sampling_method = 1;%option 1: unscented transform, option 2: random sampling with a smaller number
        if sampling_method == 1
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            n_sampling_point = 2*(obj.nx+nw);
        elseif sampling_method == 2
            n_sampling_point = 5;
            w_state = load('state_noise_small.dat');%std:0.001
            %w_state = 0.001*randn(28,62);
            %save -ascii state_noise_small.dat w_state
            w = 1;
        end
        w_avg = 1/n_sampling_point;
        K = obj.options.K;
        alpha = obj.options.alpha;
        kappa = obj.options.kappa;
         
        %initialize c and dc
        x_mean = zeros(obj.nx, obj.N);
        % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
        %c = kappa*trace(Px_init);
        %c_covariance = kappa*trace(Px_init);
        c = 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
        c_covariance = 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
        c_mean_dev = 0;
        dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
        
        % time counter
        tStart = tic;
        
        function [xdn,df] = objPlantUpdate(timestep,Sig,u_fdb_k)
            [xdot,dxdot] = obj.plant.dynamics(0,Sig,u_fdb_k);
            v_new = xdot(1:nq) + timestep*xdot(nq+1:nq+nv);
            xdot_new = [v_new;xdot(nq+1:nq+nv)];% update velocity component
            xdn = Sig + timestep*xdot_new;
            dxdot_new = full(dxdot);
            dxdot_new(1:nq,:) = [xdot(nq+1:nq+nv),zeros(nq,nx),zeros(nq,nu)] + dxdot_new(1:nq,:) + timestep*dxdot_new(nq+1:nq+nv,:);
            df = [xdot_new,eye(4),zeros(nx,nu)]+timestep*full(dxdot_new);
            % xdn = Sig + timestep*xdot;
            % df = [xdot,eye(4),zeros(4,nu)]+timestep*full(dxdot);
        end
        
        plant_update = @objPlantUpdate;
        
        noise_sample_type = 1;
        
        for k = 1:obj.N-1
            if sampling_method == 1
                %Propagate sigma points through nonlinear dynamics
                if k == 1
                    [S,d] = chol(blkdiag(Px_init, Pw),'lower');
                    if d
                        diverge = k;
                        return;
                    end
                    S = scale*S;
                    try
                        Sig(:,:,k) = [S -S];
                    catch
                        keyboard
                    end
                    
                    if isempty(w_noise)
                        for j = 1:n_sampling_point
                            Sig(1:nx,j,k) = Sig(1:nx,j,k) + x(:,k);
                        end
                    else
                        for j = 1:n_sampling_point
                            Sig(1:nx,j,k) = Sig(1:nx,j,k) + [x(:,k); w_noise(:,k)];
                        end
                    end 
                    Sig_init(:,:,k) = Sig(:,:,k);
                    %c = c + kappa*trace(Px_init);
                    %c_covariance = c_covariance + kappa*trace(Px_init);
                    %c = c + 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                    %c_covariance = c_covariance + 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                else
                    for j = 1:n_sampling_point
                        Sig_init(1:nx,j,k) = (1-alpha)*Sig(1:nx,j,k) + alpha*x(1:nx,k);
                    end
                end
            elseif sampling_method == 2
                if k == 1
                    for j = 1:n_sampling_point
                        Sig_init(:,j,k) = x(:,k) + w_state(:,j);
                    end
                    Sig(:,:,k) = Sig_init(:,:,k);
                    
                    x_mean(:,k) = zeros(nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:nx,j,k);
                    end
                    %c = c + 1/2*norm(x(:,k)-x_mean(:,k))^2;
                    %c_mean_dev = c_mean_dev + 1/2*norm(x(:,k)-x_mean(:,k))^2;
                    
                    Px_init = zeros(nx);
                    for j = 1:n_sampling_point
                        Px_init = Px_init + w*(Sig(1:nx,j,k)-x_mean(:,k))*(Sig(1:nx,j,k)-x_mean(:,k))';
                    end
                    %c = c + kappa*trace(Px_init);
                    %c_covariance = c_covariance + kappa*trace(Px_init);
                else
                    for j = 1:n_sampling_point
                        Sig_init(:,j,k) = (1-alpha)*Sig(:,j,k) + alpha*x(:,k);
                    end
                end
            end
            
            % begin of original non-parallezied version
            for j = 1:n_sampling_point
                if strcmp(obj.plant.uncertainty_source, 'physical_parameter_uncertainty')
                    % Perturb original parameter estimates with random percentage error
                    % normally distributed with standard dev = paramstd, and greater than -1
                    
                    obj.plant.m1 = 1;
                    obj.plant.m2 = 1;
                    obj.plant.l1 = 1;
                    obj.plant.l2 = 2;
                    obj.plant.b1 = 0.1;
                    obj.plant.b2 = 0.1;
                    obj.plant.lc1 = 0.5;
                    obj.plant.lc2 = 1;
                    obj.plant.Ic1 = 0.0830;
                    obj.plant.Ic2 = 0.3300;
                    
                    obj.plant.m1 = obj.plant.m1 + obj.plant.m1*param_uncertainty(j,1)/10;
                    obj.plant.m2 = obj.plant.m2 + obj.plant.m2*param_uncertainty(j,2)/10;
                    %obj.plant.l1 = obj.plant.l1 + obj.plant.l1*param_uncertainty(j,3)/10;
                    %obj.plant.l2 = obj.plant.l2 + obj.plant.l2*param_uncertainty(j,4)/10;
                    %obj.plant.b1  = obj.plant.b1 + obj.plant.b1*param_uncertainty(j,5);
                    %obj.plant.b2  = obj.plant.b2 + obj.plant.b2*param_uncertainty(j,6);
                    obj.plant.lc1 = obj.plant.lc1 + obj.plant.lc1*param_uncertainty(j,7)/10;
                    obj.plant.lc2 = obj.plant.lc2 + obj.plant.lc2*param_uncertainty(j,8)/10;
                    %obj.plant.Ic1 = obj.plant.Ic1 + obj.plant.Ic1*param_uncertainty(j,9)/10;
                    %obj.plant.Ic2 = obj.plant.Ic2 + obj.plant.Ic2*param_uncertainty(j,10)/10;
                end
                
                % add feedback control
                u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                
                if noise_sample_type == 1
                    [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k);
                elseif noise_sample_type == 2
                    xdn_analytical(:,j) = zeros(nx,1);
                    df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                    for kk=1:length(w_noise)
                        [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,0,Sig_init(1:nx,j,k),u_fdb_k);
                        xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                        df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                    end
                    xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                    df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                end
                
                % %numerical diff
                % dt = diag(max(sqrt(eps(timestep)), 1e-7));
                % dx = diag(max(sqrt(eps(Sig_init(1:nx,j,k))), 1e-7));
                % du = diag(max(sqrt(eps(u_fdb_k)),1e-7));
                %
                % [xdnp,~] = feval(plant_update,timestep+dt,Sig_init(1:nx,j,k),u_fdb_k);
                % [xdnm,~] = feval(plant_update,timestep-dt,Sig_init(1:nx,j,k),u_fdb_k);
                % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                %
                % N_finite_diff_x = length(Sig_init(1:nx,j,k));
                % for m = 1:N_finite_diff_x
                %     [xdnp,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k)+dx(:,m),u_fdb_k);
                %     [xdnm,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k)-dx(:,m),u_fdb_k);
                %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                % end
                %
                % N_finite_diff_u = length(u_fdb_k);
                % for m = 1:N_finite_diff_u
                %     [xdnp,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k+du(:,m));
                %     [xdnm,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k-du(:,m));
                %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                % end
                % %df(:,:,j) = df_numeric;
                
                % ToDo: check numerical gradient accuracy
                xdn(:,j) = xdn_analytical(:,j);
                df(:,:,j) = df_analytical(:,:,j);
                
                Sig(1:nx,j,k+1) = xdn(1:nx,j);
                dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                dfdSig(:,:,j,k+1) = (1-alpha)*df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*(1-alpha)*K;
                dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*(1-alpha)*K + alpha*df(:,2:nx+1,j);
            end
                        
            % end of original non-parallezied version
            
            % parfor jj = 1:n_sampling_point
            %     %Generate sigma points from Px(i+1)
            %     %[the sequential way to be modified]
            %     % currently, only use the initial variance matrix for the propogation
            %
            %     % a hacky way to implement the control input
            %     %[H(:,:,j),~,~,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:nxhalf,j,k),Sig(nxhalf+1:nx,j,k));
            %     %Hinv(:,:,j,k) = inv(H(:,:,j));
            %
            %     % this friction coeff samples are directly embedded
            %     % as the input argument of update() function
            %     %if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
            %     %obj.plant.uncertain_mu = w_mu(j);
            %     %end
            %
            %     % add feedback control
            %     t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
            %     %u_fdb_k(:,j) = u(:,k) - K*(Sig(1:nx,j,k) - x(:,k));
            %
            %     %[xdn,df] = obj.plant.update(t,Sig(1:nx,j,k),u_fdb_k);
            %     [xdn(:,jj),df(:,:,jj)] = feval(plant_update,timestep_updated,Sig(1:nx,jj,k),u(:,k) - K*(Sig(1:nx,jj,k) - x(:,k)));
            % end
            %
            % for jj=1:n_sampling_point
            %
            %     [H(:,:,jj),~,~,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:nxhalf,jj,k),Sig(nxhalf+1:nx,jj,k));
            %     Hinv(:,:,jj,k) = inv(H(:,:,jj));
            %
            %     Sig(1:nx,jj,k+1) = xdn(1:nx,jj);
            %
            %     dfdu(:,:,jj,k+1) = df(:,end-nu+1:end,jj);
            %     dfdSig(:,:,jj,k+1) = df(:,2:nx+1,jj) - dfdu(:,:,jj,k+1)*K;
            %     dfdx(:,:,jj,k+1) = dfdu(:,:,jj,k+1)*K;
            % end
            
            % calculate mean and variance w.r.t. [x_k] from sigma points
            x_mean(:,k+1) = zeros(nx,1);
            for j = 1:n_sampling_point
                x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
            end
            
            % % recalculate mean and variance w.r.t. [x_k] from sigma points
            % x_mean(:,k+1) = zeros(nx,1);
            % for j = 1:n_sampling_point
            %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
            % end
            
            % % check that the mean deviation term is cancelled out
            % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
            %     disp('shifting scheme is not correct')
            %     keyboard
            % end
            
            Px(:,:,k+1) = zeros(nx);
            for j = 1:n_sampling_point
                Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:nx,j,k+1)-x_mean(:,k+1))*(Sig(1:nx,j,k+1)-x_mean(:,k+1))';
            end
            
            %% ML cost
            c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*kappa*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) ...
                 + nx/2*log(2*pi);%ML mean deviation version
            %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
            c_mean_dev = c_mean_dev + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1));
            c_covariance = c_covariance + 1/2*kappa*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
            
            %% non-ML cost
            % c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
            % c = c + kappa*trace(Px(:,:,k+1));%non-ML mean deviation version
            % % note that: 1/2 coefficient is removed.
            % c_mean_dev = c_mean_dev + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;
            % c_covariance = c_covariance + kappa*trace(Px(:,:,k+1));
            
            %% non-ML version
            % % derivative of variance matrix
            % % gradient of Tr(V) w.r.t state vector x
            % dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
            % dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
            %
            % % derivative of variance matrix
            % % gradient of Tr(V) w.r.t state vector x
            % dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
            % dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
            %
            % for j=k:-1:1
            %     dTrVdx(j,:,k+1) = zeros(1,nx);
            %     dTrVdu(j,:,k+1) = zeros(1,nu);
            %
            %     % gradient w.r.t state x
            %     dSig_m_kplus1_dx_sum = zeros(nx);
            %     % gradient w.r.t control u
            %     dSig_m_kplus1_du_sum = zeros(nx,nu);
            %     dSig_i_kplus1_dx_resample = zeros(nx,nx,n_sampling_point);
            %     dSig_i_kplus1_du_resample = zeros(nx,nu,n_sampling_point);
            %
            %     for i=1:n_sampling_point
            %         if i == 1
            %             for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
            %                 % gradient of Tr(V_{k+1}) w.r.t control x and u
            %                 dSig_m_kplus1_dx = zeros(nx);
            %                 dSig_m_kplus1_du = zeros(nx,1);
            %
            %                 chain_rule_indx = k-j;
            %                 if j ~= 1
            %                     %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
            %                     dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
            %                 else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
            %                     dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2);%[proved, actually these two if-statements are same now, but derived by different ways]
            %                 end
            %                 dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
            %
            %                 while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
            %                     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
            %                     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;%note that: whether there should be a *(1-w_avg) term
            %                     chain_rule_indx = chain_rule_indx - 1;
            %                 end
            %                 dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
            %                 dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
            %             end
            %         end
            %
            %         % run 2*(nx+nw) times in total to obtain
            %         % gradient w.r.t sigma points
            %         dSig_i_kplus1_dx = zeros(nx);
            %         dSig_i_kplus1_du = zeros(nx,nu);
            %         chain_rule_indx = k-j;
            %         if j ~= 1
            %             %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
            %             dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
            %         else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
            %             dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
            %         end
            %         dSig_i_kplus1_du = dfdu(:,:,i,j+1);
            %
            %         while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
            %             dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
            %             dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
            %             chain_rule_indx = chain_rule_indx - 1;
            %         end
            %
            %         dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
            %         dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
            %
            %         %% new sigma point due to resampling mechanism
            %         %dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
            %         %dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
            %     end
            %
            %     % dSig_kplus1_dx_sum_resample = zeros(nx,nx);
            %     % dSig_kplus1_du_sum_resample = zeros(nx,nu);
            %     %
            %     % for i =1:n_sampling_point
            %     %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
            %     %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
            %     % end
            %     %
            %     % for i =1:n_sampling_point
            %     %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
            %     %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
            %     % end
            %
            %     % gradient of mean residual w.r.t state x and control u, assume norm 2
            %     dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
            %     dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
            %     %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
            %     %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
            % end
            
            %% ML version
            % derivative of variance matrix
            % gradient of Tr(V) w.r.t state vector x
            dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
            dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
            
            % derivative of variance matrix
            % gradient of Tr(V) w.r.t state vector x
            dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
            dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
            
            for j=k:-1:1
                dTrVdx(j,:,k+1) = zeros(1,nx);
                dTrVdu(j,:,k+1) = zeros(1,nu);
                
                % gradient w.r.t state x
                dSig_m_kplus1_dx_sum = zeros(nx);
                % gradient w.r.t control u
                dSig_m_kplus1_du_sum = zeros(nx,nu);
                dSig_i_kplus1_dx = zeros(nx, nx, n_sampling_point);
                dSig_i_kplus1_du = zeros(nx, nu, n_sampling_point);
                
                dSig_i_kplus1_dx_resample = zeros(nx,nx,n_sampling_point);
                dSig_i_kplus1_du_resample = zeros(nx,nu,n_sampling_point);
                
                for i=1:n_sampling_point
                    if i == 1
                        for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                            % gradient of Tr(V_{k+1}) w.r.t control x and u
                            dSig_m_kplus1_dx = zeros(nx);
                            dSig_m_kplus1_du = zeros(nx,1);
                            
                            chain_rule_indx = k-j;
                            if j ~= 1
                                %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                            else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2);%[proved, actually these two if-statements are same now, but derived by different ways]
                            end
                            dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;%note that: whether there should be a *(1-w_avg) term
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                            dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                        end
                    end
                    
                    % run 2*(nx+nw) times in total to obtain
                    % gradient w.r.t sigma points
                    %dSig_i_kplus1_dx = zeros(nx,nx,);
                    %dSig_i_kplus1_du = zeros(nx,nu);
                    chain_rule_indx = k-j;
                    if j ~= 1
                        %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                        dSig_i_kplus1_dx(:,:,i) = dfdx(:,:,i,j+1);
                    else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                        dSig_i_kplus1_dx(:,:,i) = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
                    end
                    dSig_i_kplus1_du(:,:,i) = dfdu(:,:,i,j+1);
                    
                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        dSig_i_kplus1_dx(:,:,i) = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx(:,:,i);
                        dSig_i_kplus1_du(:,:,i) = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du(:,:,i);
                        chain_rule_indx = chain_rule_indx - 1;
                    end
                    
                    %dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                    %dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                end
                
                % dSig_kplus1_dx_sum_resample = zeros(nx,nx);
                % dSig_kplus1_du_sum_resample = zeros(nx,nu);
                %
                % for i =1:n_sampling_point
                %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                % end
                %
                % for i =1:n_sampling_point
                %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                % end
                
                % gradient of mean residual w.r.t state x and control u, assume norm 2
                %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                
                dCovdx = zeros(nx,nx,nx);
                dCovdu = zeros(nx,nx,nu);
                
                for pp=1:nx
                    for i=1:n_sampling_point
                        dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                        for nn=1:nx
                            dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                            if nn == 1
                                dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                                dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                            elseif nn == nx
                                dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                                dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                            else
                                dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                                dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                                dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                            end
                            dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
                            if pp <= nu
                                dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
                            end
                        end
                    end
                    dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                    %dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                    %dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) + trace(dCovdx(:,:,pp));
                    
                    dTrVdx(j,pp,k+1) = dTrVdx(j,pp,k+1) + 1/2*kappa*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                    
                    if pp <= nu
                        dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                        %dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                        %dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) + trace(dCovdu(:,:,pp));
                        
                        dTrVdu(j,pp,k+1) = dTrVdu(j,pp,k+1) + 1/2*kappa*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                    end
                end
            end
        end

        dc = [];
        % cost gradient w.r.t x at first time step is zero
        for jj=1:obj.N % index for x_k
            dTrV_sum_dx_k = zeros(1, obj.nx);
            if (jj == 1)
                dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
            else
                %dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                dmeanR_sum_dx_k = (pinv(Px(:,:,jj)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
            end
            
            for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
            end
            dc = [dc, dmeanR_sum_dx_k+dTrV_sum_dx_k];%ML cost
            %dc = [dc, 1/2*dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];%non-ML cost
        end
        
        % cost gradient w.r.t u at first time step is zero, since
        % c(k=1) = Px(:,:,1)
        for jj=1:obj.N % index for u_k
            dTrV_sum_du_k = zeros(1, nu);
            dmeanR_sum_du_k = zeros(1, nu);
            for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
            end
            dc = [dc, dmeanR_sum_du_k+dTrV_sum_du_k];%ML cost
            %dc = [dc, 1/2*dmeanR_sum_du_k+kappa*dTrV_sum_du_k];%non-ML cost
        end
        
        % scale this robust cost
        c = obj.options.contact_robust_cost_coeff*c;
        dc = obj.options.contact_robust_cost_coeff*dc;
        c_mean_dev = obj.options.contact_robust_cost_coeff*c_mean_dev;
        c_covariance = obj.options.contact_robust_cost_coeff*c_covariance;
        fprintf('robust cost function: %4.8f\n',c);
        fprintf('robust cost mean deviation norm-2 value: %4.8f\n',c_mean_dev);
        fprintf('robust cost covariance trace value: %4.8f\n',c_covariance);
        
        for k = 1:obj.N
            q1(k) = std(Sig(1,:,k));
            q2(k) = std(Sig(2,:,k));
            v1(k) = std(Sig(3,:,k));
            v2(k) = std(Sig(4,:,k));
        end
        
        global iteration_index
        if mod(iteration_index,10) == 0
            % plot nominal model trajs
            nominal_linewidth = 1;
            color_line_type1 = 'b-';
            color_line_type2 = 'r-';
            color_line_type3 = 'k-';
            
            t = linspace(0, 6, obj.N);
            figure(10)
            hold on;
            plot(t', u, color_line_type1, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', u, color_line_type2, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('u');
            hold on;
            
            figure(22)
            subplot(2,2,1)
            hold on;
            plot(t', x(1,:), color_line_type1, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(1,:), color_line_type2, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(1,:)+q1, color_line_type3, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(1,:)-q1, color_line_type3, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('x_1')
            hold on;
            
            subplot(2,2,2)
            hold on;
            plot(t', x(2,:), color_line_type1, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(2,:), color_line_type2, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(2,:)+q2, color_line_type3, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(2,:)-q2, color_line_type3, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('x_2')
            hold on;
            
            subplot(2,2,3)
            hold on;
            plot(t', x(3,:), color_line_type1, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(3,:), color_line_type2, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(3,:)+v1, color_line_type3, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(3,:)-v1, color_line_type3, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('xdot_1')
            hold on;
            
            subplot(2,2,4)
            hold on;
            plot(t', x(4,:), color_line_type1, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(4,:), color_line_type2, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(4,:)+v2, color_line_type3, 'LineWidth',nominal_linewidth);
            hold on;
            plot(t', x_mean(4,:)-v2, color_line_type3, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('xdot_2')
            hold on;
            %keyboard
        end
        
        % figure(1)
        % kkk = 2;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_saenmpling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        %
        % figure(2)
        % kkk = 4;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_sampling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        %
        % figure(3)
        % kkk = 6;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_sampling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        %
        % figure(4)
        % kkk = 8;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_sampling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        %
        % figure(5)
        % kkk = 9;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_sampling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        %
        % figure(6)
        % kkk = 11;
        % plot(x(kkk,:),'b-');
        % hold on;
        % for j = 1:n_sampling_point
        % plot(reshape(Sig(kkk,j,:),1,[]),'r-')
        % hold on;
        % end
        
        % figure(7),hold on;plot(c_mean_dev_x,'b-');title('c_mean_dev_x');
        % figure(8),hold on;plot(c_mean_dev_xd,'b-');title('c_mean_dev_xd');
        % figure(9),hold on;plot(c_variance_x(1,:),'b-');title('c_mean_dev_x1');
        % figure(10),hold on;plot(c_variance_xd(1,:),'b-');title('c_mean_dev_xd1');
        %
        % figure(11),hold on;plot(c_mean_dev_z,'b-');title('c_mean_dev_z');
        % figure(12),hold on;plot(c_mean_dev_zd,'b-');title('c_mean_dev_zd');
        % figure(13),hold on;plot(c_variance_z(1,:),'b-');title('c_mean_dev_z1');
        % figure(14),hold on;plot(c_variance_zd(1,:),'b-');title('c_mean_dev_zd1');
        
        % figure(15)
        % clf
        % Sig_permute = permute(Sig,[1,3,2]);
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(1,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(1,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([0,7])
        % title('Sigma Point x');
        %
        % figure(16)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(7,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(7,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([0,11])
        % title('Sigma Point vx');
        %
        % figure(17)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(2,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(2,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([0,4])
        % title('Sigma Point z');
        %
        % figure(18)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(8,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(8,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-7,0])
        % title('Sigma Point vz');
        %
        % figure(19)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(3,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(3,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point first hip');
        %
        % figure(20)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(9,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(9,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point hip velocity');
        %
        % figure(21)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(4,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(4,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point first knee');
        %
        % figure(22)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(10,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(10,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point first knee velocity');
        %
        % figure(23)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(5,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(5,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point second hip');
        %
        % figure(24)
        % clf
        % for j = 1:n_sampling_point
        %     plot(Sig_permute(6,:,j),'r-')
        %     hold on;
        % end
        % hold on;
        % plot(x(6,:),'b-','Linewidth',3)
        % %xlim([0,30]);ylim([-1,1])
        % title('Sigma Point second knee');
        
        %obj.cached_Px = Px;
        %fprintf('robust cost function: %4.8f\n',c);
        
        % check gradient
        %disp('check gradient')
        c_numeric = c;
        dc_numeric = dc;
        
        X0 = [x_full; u_full];
        %X0 = X0 + randn(size(X0))*0.1;
        %
        %fun = @(X0) robustVariancecost_scaled_check(obj, X0);
        %DerivCheck(fun, X0)
        %disp('finish numerical gradient');
        
        %[c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_scaled_check(obj,X0),X0,struct('grad_method','numerical'));
        %valuecheck(dc,dc_numeric,1e-5);
        %valuecheck(c,c_numeric,1e-5);
        
        function DerivCheck(funptr, X0, ~, varargin)
            % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
            %`
            %  Checks the analytic gradient of a function 'funptr' at a point X0, and
            %  compares to numerical gradient.  Useful for checking gradients computed
            %  for fminunc and fmincon.
            %
            %  Call with same arguments as you would call for optimization (fminunc).
            %
            % $id$
            
            [~, JJ] = feval(funptr, X0, varargin{:});  % Evaluate function at X0
            
            % Pick a small vector in parameter space
            rr = sqrt(eps(X0));
            
            % Evaluate at symmetric points around X0
            f1 = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0_minus
            f2 = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0_plus
            
            % Print results
            fprintf('Derivs: Analytic vs. Finite Diff = [%.12e, %.12e]\n', dot(rr, JJ), f2-f1);
            dd = dot(rr, JJ)-f2+f1;
            fprintf('difference between numerical and analytical: %4.15f\n',dd);
        end
        
        function [c,dc] = robustVariancecost_scaled_check(obj, X0)
            x_full = X0(1:nx*obj.N);
            u_full = X0(nx*obj.N+1:end);
            
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(u_full, obj.nu, obj.N);
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.plant.getNumInputs;
            h = timestep;
            
            %compute LQR gains
            % for k = 1:obj.N-1
            %     y((k-1)*(1+nx+nu)+1+(1:nx)) = x(:,k);
            %     y((k-1)*(1+nx+nu)+1+nx+(1:nu)) = u(:,k);
            % end
            % xf = x(:,obj.N);
            % [K] = deltaLQR(obj,y,xf);
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px(:,:,1) = obj.options.Px_coeff*eye(nx);
            Px_init = Px(:,:,1);
            w_noise = [];
            Pw = [];
            
            if strcmp(obj.plant.uncertainty_source,'physical_parameter_uncertainty')
                param_uncertainty = load('physical_param_pertubation.dat');
            elseif strcmp(obj.plant.uncertainty_source,'generate_new_parameter_uncertainty_set')
                paramstd = 1/5; % Standard deviation of the parameter value percent error
                for i = 1:2*(obj.nx)
                    % Perturb original parameter estimates with random percentage error
                    % normally distributed with standard dev = paramstd, and greater than -1
                    param_uncertainty(i,:) = randn(1,10)*paramstd;
                    while sum(param_uncertainty(i,:)<=-1)~=0
                        param_uncertainty(param_uncertainty(i,:)<-1) = randn(1,sum(param_uncertainty(i,:)<-1))*paramstd;
                    end
                end
                save -ascii physical_param_pertubation.dat param_uncertainty
                w_noise = param_uncertainty;
            end
            
            % disturbance variance
            % currently only consider object horizontal 2D position and friction coefficient
            nw = 0;%size(Pw,1);
            sampling_method = 1;%option 1: unscented transform, option 2: random sampling with a smaller number
            if sampling_method == 1
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
                n_sampling_point = 2*(obj.nx+nw);
            elseif sampling_method == 2
                n_sampling_point = 5;
                w_state = load('state_noise_small.dat');%std:0.001
                %w_state = 0.001*randn(28,62);
                %save -ascii state_noise_small.dat w_state
                w = 1;
            end
            w_avg = 1/n_sampling_point;
            K = obj.options.K;
            alpha = obj.options.alpha;
            kappa = obj.options.kappa;
            
            %initialize c and dc
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            %c = 0;
            %c = kappa*trace(Px_init);
            c = 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
            %c_covariance = kappa*trace(Px_init);
            c_covariance = 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
            c_mean_dev = 0;
            dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(timestep,Sig,u_fdb_k)
                [xdot,dxdot] = obj.plant.dynamics(0,Sig,u_fdb_k);
                v_new = xdot(1:nq) + timestep*xdot(nq+1:nq+nv);
                xdot_new = [v_new;xdot(nq+1:nq+nv)];% update velocity component
                xdn = Sig + timestep*xdot_new;
                dxdot_new = full(dxdot);
                dxdot_new(1:nq,:) = [xdot(nq+1:nq+nv),zeros(nq,nx),zeros(nq,nu)] + dxdot_new(1:nq,:) + timestep*dxdot_new(nq+1:nq+nv,:);
                df = [xdot_new,eye(4),zeros(nx,nu)]+timestep*full(dxdot_new);
            end
            
            plant_update = @objPlantUpdate;
            
            noise_sample_type = 1;
            
            for k = 1:obj.N-1
                if sampling_method == 1
                    %Propagate sigma points through nonlinear dynamics
                    if k == 1
                        [S,d] = chol(blkdiag(Px_init, Pw),'lower');
                        if d
                            diverge = k;
                            return;
                        end
                        S = scale*S;
                        try
                            Sig(:,:,k) = [S -S];
                        catch
                            keyboard
                        end
                        
                        if isempty(w_noise)
                            for j = 1:n_sampling_point
                                Sig(1:nx,j,k) = Sig(1:nx,j,k) + x(:,k);
                            end
                        else
                            for j = 1:n_sampling_point
                                Sig(1:nx,j,k) = Sig(1:nx,j,k) + [x(:,k); w_noise(:,k)];
                            end
                        end
                        Sig_init(:,:,k) = Sig(:,:,k);
                        %c = c + kappa*trace(Px_init);
                        %c_covariance = c_covariance + kappa*trace(Px_init);
                        c = c + 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                        c_covariance = c_covariance + 1/2*kappa*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                    else
                        for j = 1:n_sampling_point
                            Sig_init(1:nx,j,k) = (1-alpha)*Sig(1:nx,j,k) + alpha*x(1:nx,k);
                        end
                    end
                elseif sampling_method == 2
                    if k == 1
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = x(:,k) + w_state(:,j);
                        end
                        Sig(:,:,k) = Sig_init(:,:,k);
                        
                        x_mean(:,k) = zeros(nx,1);
                        for j = 1:n_sampling_point
                            x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:nx,j,k);
                        end
                        c = c + 1/2*norm(x(:,k)-x_mean(:,k))^2;
                        c_mean_dev = c_mean_dev + 1/2*norm(x(:,k)-x_mean(:,k))^2;
                        
                        Px_init = zeros(nx);
                        for j = 1:n_sampling_point
                            Px_init = Px_init + w*(Sig(1:nx,j,k)-x_mean(:,k))*(Sig(1:nx,j,k)-x_mean(:,k))';
                        end
                        c = c + kappa*trace(Px_init);
                        c_covariance = c_covariance + kappa*trace(Px_init);
                    else
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = (1-alpha)*Sig(:,j,k) + alpha*x(:,k);
                        end
                    end
                end
                
                % begin of original non-parallezied version
                for j = 1:n_sampling_point
                    if strcmp(obj.plant.uncertainty_source, 'physical_parameter_uncertainty')
                        % Perturb original parameter estimates with random percentage error
                        % normally distributed with standard dev = paramstd, and greater than -1
                        
                        obj.plant.m1 = 1;
                        obj.plant.m2 = 1;
                        obj.plant.l1 = 1;
                        obj.plant.l2 = 2;
                        obj.plant.b1 = 0.1;
                        obj.plant.b2 = 0.1;
                        obj.plant.lc1 = 0.5;
                        obj.plant.lc2 = 1;
                        obj.plant.Ic1 = 0.0830;
                        obj.plant.Ic2 = 0.3300;
                        
                        obj.plant.m1 = obj.plant.m1 + obj.plant.m1*param_uncertainty(j,1)/10;
                        obj.plant.m2 = obj.plant.m2 + obj.plant.m2*param_uncertainty(j,2)/10;
                        %obj.plant.l1 = obj.plant.l1 + obj.plant.l1*param_uncertainty(j,3)/10;
                        %obj.plant.l2 = obj.plant.l2 + obj.plant.l2*param_uncertainty(j,4)/10;
                        %obj.plant.b1  = obj.plant.b1 + obj.plant.b1*param_uncertainty(j,5);
                        %obj.plant.b2  = obj.plant.b2 + obj.plant.b2*param_uncertainty(j,6);
                        obj.plant.lc1 = obj.plant.lc1 + obj.plant.lc1*param_uncertainty(j,7)/10;
                        obj.plant.lc2 = obj.plant.lc2 + obj.plant.lc2*param_uncertainty(j,8)/10;
                        %obj.plant.Ic1 = obj.plant.Ic1 + obj.plant.Ic1*param_uncertainty(j,9)/10;
                        %obj.plant.Ic2 = obj.plant.Ic2 + obj.plant.Ic2*param_uncertainty(j,10)/10;
                    end
                    
                    % add feedback control
                    u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                    
                    if noise_sample_type == 1
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k);
                    elseif noise_sample_type == 2
                        xdn_analytical(:,j) = zeros(nx,1);
                        df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                        for kk=1:length(w_noise)
                            [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,0,Sig_init(1:nx,j,k),u_fdb_k);
                            xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                            df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                        end
                        xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                        df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                    end
                    
                    % %numerical diff
                    % dt = diag(max(sqrt(eps(timestep)), 1e-7));
                    % dx = diag(max(sqrt(eps(Sig_init(1:nx,j,k))), 1e-7));
                    % du = diag(max(sqrt(eps(u_fdb_k)),1e-7));
                    %
                    % [xdnp,~] = feval(plant_update,timestep+dt,Sig_init(1:nx,j,k),u_fdb_k);
                    % [xdnm,~] = feval(plant_update,timestep-dt,Sig_init(1:nx,j,k),u_fdb_k);
                    % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                    %
                    % N_finite_diff_x = length(Sig_init(1:nx,j,k));
                    % for m = 1:N_finite_diff_x
                    %     [xdnp,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k)+dx(:,m),u_fdb_k);
                    %     [xdnm,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k)-dx(:,m),u_fdb_k);
                    %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                    % end
                    %
                    % N_finite_diff_u = length(u_fdb_k);
                    % for m = 1:N_finite_diff_u
                    %     [xdnp,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k+du(:,m));
                    %     [xdnm,~] = feval(plant_update,timestep,Sig_init(1:nx,j,k),u_fdb_k-du(:,m));
                    %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                    % end
                    % %df(:,:,j) = df_numeric;
                    
                    % ToDo: check numerical gradient accuracy
                    xdn(:,j) = xdn_analytical(:,j);
                    df(:,:,j) = df_analytical(:,:,j);
                    
                    Sig(1:nx,j,k+1) = xdn(1:nx,j);
                    dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                    dfdSig(:,:,j,k+1) = (1-alpha)*df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*(1-alpha)*K;
                    dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*(1-alpha)*K + alpha*df(:,2:nx+1,j);
                end
                
                % end of original non-parallezied version
                
                % parfor jj = 1:n_sampling_point
                %     %Generate sigma points from Px(i+1)
                %     %[the sequential way to be modified]
                %     % currently, only use the initial variance matrix for the propogation
                %
                %     % a hacky way to implement the control input
                %     %[H(:,:,j),~,~,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:nxhalf,j,k),Sig(nxhalf+1:nx,j,k));
                %     %Hinv(:,:,j,k) = inv(H(:,:,j));
                %
                %     % this friction coeff samples are directly embedded
                %     % as the input argument of update() function
                %     %if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                %     %obj.plant.uncertain_mu = w_mu(j);
                %     %end
                %
                %     % add feedback control
                %     t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                %     %u_fdb_k(:,j) = u(:,k) - K*(Sig(1:nx,j,k) - x(:,k));
                %
                %     %[xdn,df] = obj.plant.update(t,Sig(1:nx,j,k),u_fdb_k);
                %     [xdn(:,jj),df(:,:,jj)] = feval(plant_update,timestep_updated,Sig(1:nx,jj,k),u(:,k) - K*(Sig(1:nx,jj,k) - x(:,k)));
                % end
                %
                % for jj=1:n_sampling_point
                %
                %     [H(:,:,jj),~,~,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:nxhalf,jj,k),Sig(nxhalf+1:nx,jj,k));
                %     Hinv(:,:,jj,k) = inv(H(:,:,jj));
                %
                %     Sig(1:nx,jj,k+1) = xdn(1:nx,jj);
                %
                %     dfdu(:,:,jj,k+1) = df(:,end-nu+1:end,jj);
                %     dfdSig(:,:,jj,k+1) = df(:,2:nx+1,jj) - dfdu(:,:,jj,k+1)*K;
                %     dfdx(:,:,jj,k+1) = dfdu(:,:,jj,k+1)*K;
                % end
                
                % calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                end
                
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                % end
                
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                
                Px(:,:,k+1) = zeros(nx);
                for j = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:nx,j,k+1)-x_mean(:,k+1))*(Sig(1:nx,j,k+1)-x_mean(:,k+1))';
                end
                
                % accumulate returned cost
                c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*kappa*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) ...
                    + nx/2*log(2*pi);%ML mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                c_mean_dev = c_mean_dev + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1));
                c_covariance = c_covariance + 1/2*kappa*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                
                %c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                %c = c + kappa*trace(Px(:,:,k+1));%non-ML mean deviation version
                % note that: 1/2 coefficient is removed.
                %c_mean_dev = c_mean_dev + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;
                %c_covariance = c_covariance + kappa*trace(Px(:,:,k+1));
                
                %             %% non-ML version
                %             % derivative of variance matrix
                %             % gradient of Tr(V) w.r.t state vector x
                %             dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                %             dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                %
                %             % derivative of variance matrix
                %             % gradient of Tr(V) w.r.t state vector x
                %             dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
                %             dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                %
                %             for j=k:-1:1
                %                 dTrVdx(j,:,k+1) = zeros(1,nx);
                %                 dTrVdu(j,:,k+1) = zeros(1,nu);
                %
                %                 % gradient w.r.t state x
                %                 dSig_m_kplus1_dx_sum = zeros(nx);
                %                 % gradient w.r.t control u
                %                 dSig_m_kplus1_du_sum = zeros(nx,nu);
                %                 dSig_i_kplus1_dx_resample = zeros(nx,nx,n_sampling_point);
                %                 dSig_i_kplus1_du_resample = zeros(nx,nu,n_sampling_point);
                %
                %                 for i=1:n_sampling_point
                %                     if i == 1
                %                         for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                %                             % gradient of Tr(V_{k+1}) w.r.t control x and u
                %                             dSig_m_kplus1_dx = zeros(nx);
                %                             dSig_m_kplus1_du = zeros(nx,1);
                %
                %                             chain_rule_indx = k-j;
                %                             if j ~= 1
                %                                 %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                %                                 dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                %                             else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                %                                 dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2);%[proved, actually these two if-statements are same now, but derived by different ways]
                %                             end
                %                             dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                %
                %                             while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                %                                 dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                %                                 dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;%note that: whether there should be a *(1-w_avg) term
                %                                 chain_rule_indx = chain_rule_indx - 1;
                %                             end
                %                             dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                %                             dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                %                         end
                %                     end
                %
                %                     % run 2*(nx+nw) times in total to obtain
                %                     % gradient w.r.t sigma points
                %                     dSig_i_kplus1_dx = zeros(nx);
                %                     dSig_i_kplus1_du = zeros(nx,nu);
                %                     chain_rule_indx = k-j;
                %                     if j ~= 1
                %                         %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                %                         dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                %                     else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                %                         dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
                %                     end
                %                     dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                %
                %                     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                %                         dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                %                         dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                %                         chain_rule_indx = chain_rule_indx - 1;
                %                     end
                %
                %                     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                %                     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                %
                %                     %% new sigma point due to resampling mechanism
                %                     %dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                %                     %dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                %                 end
                %
                %                 % dSig_kplus1_dx_sum_resample = zeros(nx,nx);
                %                 % dSig_kplus1_du_sum_resample = zeros(nx,nu);
                %                 %
                %                 % for i =1:n_sampling_point
                %                 %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                %                 %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                %                 % end
                %                 %
                %                 % for i =1:n_sampling_point
                %                 %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                %                 %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                %                 % end
                %
                %                 % gradient of mean residual w.r.t state x and control u, assume norm 2
                %                 dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                %                 dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                %                 %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                %                 %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                %             end
                
                %% ML version
                % derivative of variance matrix
                % gradient of Tr(V) w.r.t state vector x
                dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                % derivative of variance matrix
                % gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
                dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                for j=k:-1:1
                    dTrVdx(j,:,k+1) = zeros(1,nx);
                    dTrVdu(j,:,k+1) = zeros(1,nu);
                    
                    % gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(nx);
                    % gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(nx,nu);
                    dSig_i_kplus1_dx = zeros(nx, nx, n_sampling_point);
                    dSig_i_kplus1_du = zeros(nx, nu, n_sampling_point);
                    
                    dSig_i_kplus1_dx_resample = zeros(nx,nx,n_sampling_point);
                    dSig_i_kplus1_du_resample = zeros(nx,nu,n_sampling_point);
                    
                    for i=1:n_sampling_point
                        if i == 1
                            for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                % gradient of Tr(V_{k+1}) w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(nx);
                                dSig_m_kplus1_du = zeros(nx,1);
                                
                                chain_rule_indx = k-j;
                                if j ~= 1
                                    %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2);%[proved, actually these two if-statements are same now, but derived by different ways]
                                end
                                dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                
                                while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;%note that: whether there should be a *(1-w_avg) term
                                    chain_rule_indx = chain_rule_indx - 1;
                                end
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run 2*(nx+nw) times in total to obtain
                        % gradient w.r.t sigma points
                        %dSig_i_kplus1_dx = zeros(nx,nx,);
                        %dSig_i_kplus1_du = zeros(nx,nu);
                        chain_rule_indx = k-j;
                        if j ~= 1
                            %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                            dSig_i_kplus1_dx(:,:,i) = dfdx(:,:,i,j+1);
                        else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            dSig_i_kplus1_dx(:,:,i) = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
                        end
                        dSig_i_kplus1_du(:,:,i) = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx(:,:,i) = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx(:,:,i);
                            dSig_i_kplus1_du(:,:,i) = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du(:,:,i);
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                        %dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                        %dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        
                        %% new sigma point due to resampling mechanism
                        %dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                        %dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                    end
                    
                    % dSig_kplus1_dx_sum_resample = zeros(nx,nx);
                    % dSig_kplus1_du_sum_resample = zeros(nx,nu);
                    %
                    % for i =1:n_sampling_point
                    %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                    %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                    % end
                    %
                    % for i =1:n_sampling_point
                    %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                    %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                    % end
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                    %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    
                    dCovdx = zeros(nx,nx,nx);
                    dCovdu = zeros(nx,nx,nu);
                    
                    for pp=1:nx
                        for i=1:n_sampling_point
                            dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                            for nn=1:nx
                                dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                                if nn == 1
                                    dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                                    dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                                elseif nn == nx
                                    dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                                    dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                                else
                                    dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                                    dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                                    dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                                end
                                dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
                                if pp <= nu
                                    dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
                                end
                            end
                        end
                        dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                        %dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                        %dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) + trace(dCovdx(:,:,pp));
                        
                        dTrVdx(j,pp,k+1) = dTrVdx(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdx(:,:,pp));
                        
                        if pp <= nu
                            dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                            %dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                            %dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) + trace(dCovdu(:,:,pp));
                            
                            dTrVdu(j,pp,k+1) = dTrVdu(j,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*dCovdu(:,:,pp));
                        end
                    end
                end
                
                % induced by resampling mechanism
                % dTrVdx(k+1,:,k+1) = zeros(nx); since
                % dSig_i_kplus1_dx_resample(:,:,k+1) = zeros(nx);
            end
            
            dc = [];
            % cost gradient w.r.t x at first time step is zero
            for jj=1:obj.N % index for x_k
                dTrV_sum_dx_k = zeros(1, obj.nx);
                if (jj == 1)
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                else
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                    %dmeanR_sum_dx_k = (inv(Px(:,:,jj))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                end
                
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                    dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
                end
                dc = [dc, 1/2*dmeanR_sum_dx_k+dTrV_sum_dx_k];
                %dc = [dc, dmeanR_sum_dx_k];
                %dc = [dc, dTrV_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            % c(k=1) = Px(:,:,1)
            for jj=1:obj.N % index for u_k
                dTrV_sum_du_k = zeros(1, nu);
                dmeanR_sum_du_k = zeros(1, nu);
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                end
                dc = [dc, 1/2*dmeanR_sum_du_k+dTrV_sum_du_k];
                %dc = [dc, dmeanR_sum_du_k];
                %dc = [dc, dTrV_sum_du_k];
            end
            
            % scale this robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
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
    
    % discrete-time LQR modified from Zac's Dirtrel code
    function [K] = deltaLQR(obj,y,xf)
        nX = obj.nX;
        nU = obj.nU;
        nW = obj.nW;
        N = obj.N;
        
        %check to see if we actually need to do anything
        
        %Get dynamics derivatives along trajectory
        A = zeros(nX,nX,N-1);
        B = zeros(nX,nU,N-1);
        %G = zeros(nX,nW,N-1);
        
        dA = zeros(nX*nX,2*(1+nX+nU),N-1);
        dB = zeros(nX*nU,2*(1+nX+nU),N-1);
        %dG = zeros(nX*nW,2*(1+nX+nU),N-1);
        for k = 1:(N-2)
            [~,dx,d2x] = obj.dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+y((k)*(1+nX+nU)+1+(1:nX))),y((k-1)*(1+nX+nU)+1+nX+(1:nU)));
            A(:,:,k) = dx(:,1+(1:nX));
            B(:,:,k) = dx(:,1+nX+(1:nU));
            %G(:,:,k) = dx(:,1+nX+nU+(1:nW));
            dvec = reshape(d2x,nX*(1+nX+nU),1+nX+nU);
            dA(:,:,k) = [dvec(nX+(1:nX*nX),1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), dvec(nX+(1:nX*nX),1+nX+(1:nU)), zeros(nX*nX,1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), zeros(nX*nX,nU)];
            dB(:,:,k) = [dvec((1+nX)*nX+(1:nX*nU),1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), dvec((1+nX)*nX+(1:nX*nU),1+nX+(1:nU)), zeros(nX*nU,1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), zeros(nX*nU,nU)];
            %dG(:,:,k) = [dvec((1+nX+nU)*nX+(1:nX*nW),1), .5*dvec((1+nX+nU)*nX+(1:nX*nW),1+(1:nX)), dvec((1+nX+nU)*nX+(1:nX*nW),1+nX+(1:nU)), zeros(nX*nW,1), .5*dvec((1+nX+nU)*nX+(1:nX*nW),1+(1:nX)), zeros(nX*nW,nU)];
        end
        k = N-1;
        [~,dx,d2x] = obj.dynamics(y((k-1)*(1+nX+nU)+1),.5*(y((k-1)*(1+nX+nU)+1+(1:nX))+xf'),y((k-1)*(1+nX+nU)+1+nX+(1:nU)));
        A(:,:,k) = dx(:,1+(1:nX));
        B(:,:,k) = dx(:,1+nX+(1:nU));
        %G(:,:,k) = dx(:,1+nX+nU+(1:nW));
        dvec = reshape(d2x,nX*(1+nX+nU),1+nX+nU);
        dA(:,:,k) = [dvec(nX+(1:nX*nX),1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), dvec(nX+(1:nX*nX),1+nX+(1:nU)), zeros(nX*nX,1), .5*dvec(nX+(1:nX*nX),1+(1:nX)), zeros(nX*nX,nU)];
        dB(:,:,k) = [dvec((1+nX)*nX+(1:nX*nU),1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), dvec((1+nX)*nX+(1:nX*nU),1+nX+(1:nU)), zeros(nX*nU,1), .5*dvec((1+nX)*nX+(1:nX*nU),1+(1:nX)), zeros(nX*nU,nU)];
        %dG(:,:,k) = [dvec((1+nX+nU)*nX+(1:nX*nW),1), .5*dvec((1+nX+nU)*nX+(1:nX*nW),1+(1:nX)), dvec((1+nX+nU)*nX+(1:nX*nW),1+nX+(1:nU)), zeros(nX*nW,1), .5*dvec((1+nX+nU)*nX+(1:nX*nW),1+(1:nX)), zeros(nX*nW,nU)];
        
        %Solve Riccati Equation
        P = 100*eye(4);
        obj.Q = diag([10 10 1 1]);
        obj.R = .1;
        
        dP = zeros(nX*nX,(N-1)*(1+nX+nU)+nX);
        K = zeros(nU,nX,N-1);
        dK = zeros(nU*nX,(N-1)*(1+nX+nU)+nX,N-1);
        
        k = N-1;
        K(:,:,k) = (B(:,:,k).'*P*B(:,:,k)+obj.R)\(B(:,:,k).'*P*A(:,:,k));
        dKdA = kron(eye(nX),(B(:,:,k)'*P*B(:,:,k)+obj.R)\B(:,:,k)'*P);
        dKdB = kron(A(:,:,k)'*P, inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*comm(nX,nU) - kron(A(:,:,k)'*P*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*P*B(:,:,k)+obj.R)', inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*(kron(eye(nU), B(:,:,k)'*P) + kron(B(:,:,k)'*P, eye(nU))*comm(nX,nU));
        dKdP = kron(A(:,:,k)', (B(:,:,k)'*P*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*P*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*P*B(:,:,k)+obj.R)', inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
        dK(:,:,k) = dKdP*dP;
        dK(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU)),k) = dKdA*dA(:,1:(1+nX+nU),k) + dKdB*dB(:,1:(1+nX+nU),k);
        dK(:,k*(1+nX+nU)+(1:nX),k) = dKdA*dA(:,(1+nX+nU+1)+(1:nX),k) + dKdB*dB(:,(1+nX+nU+1)+(1:nX),k);
        dPdA = kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*comm(nX,nX);
        dPdB = -kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P)*kron(K(:,:,k)', eye(nX)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*kron(eye(nX), K(:,:,k)')*comm(nX,nU);
        dPdK = kron(eye(nX), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nX))*comm(nU,nX) - kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P)*kron(eye(nX), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*kron(B(:,:,k), eye(nX))*comm(nU,nX);
        dPdP = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
        P = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*P*(A(:,:,k) - B(:,:,k)*K(:,:,k));
        dP = dPdP*dP + dPdK*dK(:,:,k);
        dP(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) = dP(:,(k-1)*(1+nX+nU)+(1:(1+nX+nU))) + dPdA*dA(:,1:(1+nX+nU),k) + dPdB*dB(:,1:(1+nX+nU),k);
        dP(:,k*(1+nX+nU)+(1:nX)) = dP(:,k*(1+nX+nU)+(1:nX))+ dPdA*dA(:,(1+nX+nU+1)+(1:nX),k) + dPdB*dB(:,(1+nX+nU+1)+(1:nX),k);
        for k = (N-2):-1:1
            K(:,:,k) = (B(:,:,k).'*P*B(:,:,k)+obj.R)\(B(:,:,k).'*P*A(:,:,k));
            dKdA = kron(eye(nX),(B(:,:,k)'*P*B(:,:,k)+obj.R)\B(:,:,k)'*P);
            dKdB = kron(A(:,:,k)'*P, inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*comm(nX,nU) - kron(A(:,:,k)'*P*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*P*B(:,:,k)+obj.R)', inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*(kron(eye(nU), B(:,:,k)'*P) + kron(B(:,:,k)'*P, eye(nU))*comm(nX,nU));
            dKdP = kron(A(:,:,k)', (B(:,:,k)'*P*B(:,:,k)+obj.R)\B(:,:,k)') - kron(A(:,:,k)'*P*B(:,:,k), eye(nU))*kron(inv(B(:,:,k)'*P*B(:,:,k)+obj.R)', inv(B(:,:,k)'*P*B(:,:,k)+obj.R))*kron(B(:,:,k)', B(:,:,k)');
            dK(:,:,k) = dKdP*dP;
            dK(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU)),k) = dK(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU)),k) + dKdA*dA(:,:,k) + dKdB*dB(:,:,k);
            
            dPdA = kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P) + kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*comm(nX,nX);
            dPdB = -kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P)*kron(K(:,:,k)', eye(nX)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*kron(eye(nX), K(:,:,k)')*comm(nX,nU);
            dPdK = kron(eye(nX), K(:,:,k)'*obj.R) + kron(K(:,:,k)'*obj.R, eye(nX))*comm(nU,nX) - kron(eye(nX), (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P)*kron(eye(nX), B(:,:,k)) - kron((A(:,:,k)-B(:,:,k)*K(:,:,k))'*P, eye(nX))*kron(B(:,:,k), eye(nX))*comm(nU,nX);
            dPdP = kron((A(:,:,k)-B(:,:,k)*K(:,:,k))', (A(:,:,k)-B(:,:,k)*K(:,:,k))');
            
            P = obj.Q + K(:,:,k).'*obj.R*K(:,:,k) + (A(:,:,k) - B(:,:,k)*K(:,:,k)).'*P*(A(:,:,k) - B(:,:,k)*K(:,:,k));
            dP = dPdP*dP + dPdK*dK(:,:,k);
            dP(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) = dP(:,(k-1)*(1+nX+nU)+(1:2*(1+nX+nU))) + dPdA*dA(:,:,k) + dPdB*dB(:,:,k);
        end
        
%         obj.K_handle.data = K;
        return
    end
    
    function [f,df,d2f] = robust_dynamics(obj,h,x,u,w)
      % Euler integration of continuous dynamics
      if nargout == 1
          xdot = obj.plant.dynamics_w(0,x,u,w);
          f = x + h*xdot;
      elseif nargout == 2
        [xdot,dxdot] = obj.plant.dynamics_w(0,x,u,w);
        f = x + h*xdot;
        df = [xdot ... h
          eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
          h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
          h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
      else %nargout == 3
          [xdot,dxdot,d2xdot] = obj.plant.dynamics_w(0,x,u,w);
          f = x + h*xdot;
          df = [xdot ... h
            eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
            h*dxdot(:,1+obj.nX+(1:obj.nU)) ... u
            h*dxdot(:,1+obj.nX+obj.nU+(1:obj.nW))]; % w
          d2f = h*d2xdot;
          d2f(:,1:(1+obj.nX+obj.nU+obj.nW)) = dxdot;
          d2f(:,1:(1+obj.nX+obj.nU+obj.nW):end) = dxdot;
      end
    end
    
    function [f,df,d2f] = dynamics(obj,h,x,u)
        x = x';
        u = u';
        % Euler integration of continuous dynamics
        if nargout == 1
            xdot = obj.plant.dynamics(0,x,u);
            f = x + h*xdot;
        elseif nargout == 2
            [xdot,dxdot] = obj.plant.dynamics(0,x,u);
            f = x + h*xdot;
            df = [xdot ... h
                eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
                h*dxdot(:,1+obj.nX+(1:obj.nU))]; % u
        else %nargout == 3
            [xdot,dxdot,d2xdot] = obj.plant.dynamics(0,x,u);
            f = x + h*xdot;
            df = [xdot ... h
                eye(obj.nX) + h*dxdot(:,1+(1:obj.nX)) ... x0
                h*dxdot(:,1+obj.nX+(1:obj.nU))]; % u
            d2f = h*d2xdot;
            d2f(:,1:(1+obj.nX+obj.nU)) = dxdot;
            d2f(:,1:(1+obj.nX+obj.nU):end) = dxdot;
        end
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