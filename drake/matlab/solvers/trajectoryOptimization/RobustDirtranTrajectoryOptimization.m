classdef RobustDirtranTrajectoryOptimization < DirectTrajectoryOptimization
    % Direct transcription trajectory optimization
    %  implements multiple possible integration schemes for the dynamics
    %  constraints xdot = f(x,u) and for for integrating the running cost
    %  term.
    %
    %  For forward euler integratino:
    %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(x(k),u(k))
    %    integrated cost is sum of g(h(k),x(k),u(k))
    %  For backward euler integration:
    %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(x(k+1),u(k))
    %    integrated cost is sum of g(h(k),x(k+1),u(k))
    %  For midpoint integration:
    %    dynamics constraints are: x(k+1) = x(k) + h(k)*f(.5*x(k)+.5*x(k+1),.5*u(k)+.5*u(k+1))
    %    integrated cost is sum of g(h(k),.5*x(k)+.5*x(k+1),.5*u(k)+.5*u(k+1))
    properties (Constant)
        FORWARD_EULER = 1;
        BACKWARD_EULER = 2;
        MIDPOINT = 3;  % DEFAULT
    end
    
    properties
        nq
        nv
        nu
        nx
        h
    end
    
    methods
        function obj = RobustDirtranTrajectoryOptimization(plant,N,duration,options)
            if nargin < 4
                options = struct();
            end
            if ~isfield(options,'integration_method')
                options.integration_method = RobustDirtranTrajectoryOptimization.MIDPOINT;
            end
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
        end
        
        function obj = addDynamicConstraints(obj)
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            N = obj.N;
            obj.nx = nX;
            obj.nq = nX/2;
            obj.nv = nX/2;
            obj.nu = nU;
            
            constraints = cell(N-1,1);
            dyn_inds = cell(N-1,1);
            
            switch obj.options.integration_method
                case RobustDirtranTrajectoryOptimization.FORWARD_EULER
                    n_vars = 2*nX + nU + 1;
                    cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.forward_constraint_fun);
                case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
                    n_vars = 2*nX + nU + 1;
                    cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.backward_constraint_fun);
                case RobustDirtranTrajectoryOptimization.MIDPOINT
                    n_vars = 2*nX + 2*nU + 1;
                    cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.midpoint_constraint_fun);
                otherwise
                    error('Drake:RobustDirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
            end
            
            for i=1:obj.N-1,
                switch obj.options.integration_method
                    case RobustDirtranTrajectoryOptimization.FORWARD_EULER
                        dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
                        dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    case RobustDirtranTrajectoryOptimization.MIDPOINT
                        dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
                    otherwise
                        error('Drake:RobustDirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
                end
                constraints{i} = cnstr;
                
                obj = obj.addConstraint(constraints{i}, dyn_inds{i});
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
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px(:,:,1) = obj.options.Px_coeff*eye(nx);
            Px_init = Px(:,:,1);
            w_noise = [];
            Pw = [];
            
%             if strcmp(obj.plant.uncertainty_source,'friction_coeff')
%                 w_mu = obj.plant.uncertain_mu_set;
%                 w_noise = [w_mu];
%                 Pw = diag([0.01]);
%             elseif strcmp(obj.plant.uncertainty_source,'object_initial_position')
%                 w_phi = obj.plant.uncertain_position_set;
%                 w_noise = [w_phi];
%                 Pw = diag([0.0032,0.0037]);
%             elseif strcmp(obj.plant.uncertainty_source,'object_initial_orientation')
%                 w_ori = obj.plant.uncertain_orientation_set;
%                 w_noise = [w_ori];
%                 Pw = diag([0.025,0.025]);
%             elseif strcmp(obj.plant.uncertainty_source,'friction_coeff+object_initial_position')
%                 w_mu = obj.plant.uncertain_mu_set;
%                 w_phi = obj.plant.uncertain_position_set;
%                 w_noise = [w_mu;w_phi];
%                 Pw = diag([0.01, 0.0032,0.0037]);
%             elseif strcmp(obj.plant.uncertainty_source,'friction_coeff+object_initial_orientation')
%                 w_mu = obj.plant.uncertain_mu_set;
%                 w_ori = obj.plant.uncertain_orientation_set;
%                 w_noise = [w_mu;w_ori];
%                 Pw = diag([0.01, 0.025,0.025]);
%             elseif isempty(obj.plant.uncertainty_source)
%                 Pw = [];
%                 w_noise = [];
%             elseif strcmp(obj.plant.uncertainty_source,'generate_new_noise_set')
%                 w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
%                 save -ascii friction_coeff_noise.dat w_mu
%                 %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
%                 %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(3,3));
%                 %w_phi = [x;y];%height noise
%                 %save -ascii initial_position_noise.dat w_phi
%                 w_noise = [w_mu];
%             end
            
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
            c = trace(Px_init);
            %c = 1/2*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
            c_covariance = 0;
            c_mean_dev = 0;
            dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(timestep,Sig,u_fdb_k)
                [xdot,dxdot] = obj.plant.dynamics(0,Sig,u_fdb_k);
                xdn = Sig + timestep*xdot;
                df = [xdot,eye(4),zeros(4,nu)]+timestep*full(dxdot);
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
                        c = c + kappa*trace(Px_init);
                        c_covariance = c_covariance + kappa*trace(Px_init);
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
%                     if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
%                         obj.plant.uncertain_mu = w_mu(j);
%                     elseif strcmp(obj.plant.uncertainty_source, 'object_initial_position')
%                         obj.plant.uncertain_phi = w_phi(:,j);
%                         Sig_init(9:10,j,k) = Sig_init(9:10,j,k) + obj.plant.uncertain_phi;%object x and y position uncertainty
%                     elseif strcmp(obj.plant.uncertainty_source, 'object_initial_orientation')
%                         obj.plant.uncertain_ori = w_ori(:,j);
%                         Sig_init(12:13,j,k) = Sig_init(12:13,j,k) + obj.plant.uncertain_ori;%object x and y orientation uncertainty
%                     elseif strcmp(obj.plant.uncertainty_source, 'friction_coeff+object_initial_position')
%                         obj.plant.uncertain_mu = w_mu(j);
%                         obj.plant.uncertain_phi = w_phi(:,j);
%                         Sig_init(9:10,j,k) = Sig_init(9:10,j,k) + obj.plant.uncertain_phi;%object x and y position uncertainty
%                     elseif strcmp(obj.plant.uncertainty_source, 'friction_coeff+object_initial_orientation')
%                         obj.plant.uncertain_mu = w_mu(j);
%                         obj.plant.uncertain_ori = w_ori(:,j);
%                         Sig_init(12:13,j,k) = Sig_init(12:13,j,k) + obj.plant.uncertain_ori;%object x and y orientation uncertainty
%                     end
                    
                    %[H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:nx/2,j,k),Sig_init(nx/2+1:nx,j,k));
                    %Hinv(:,:,j,k) = inv(H);
                    
                    % add feedback control
                    u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                    
                    external_force_index = j;
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
                    %df(:,:,j) = df_numeric;
                    
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
                c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) ...
                %    + nx/2*log(2*pi);%ML mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                %c = c + (x(:,k+1)-x_mean(:,k+1))'*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                c = c + kappa*trace(Px(:,:,k+1));%ML mean deviation version
                % note that: 1/2 coefficient is removed.
                c_mean_dev = c_mean_dev + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;
                c_covariance = c_covariance + kappa*trace(Px(:,:,k+1));
                
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
                        dSig_i_kplus1_dx = zeros(nx);
                        dSig_i_kplus1_du = zeros(nx,nu);
                        chain_rule_indx = k-j;
                        if j ~= 1
                            %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                        else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
                        end
                        dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                            dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                        dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                        dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        
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
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                    %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
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
            c_mean_dev = obj.options.contact_robust_cost_coeff*c_mean_dev;
            c_covariance = obj.options.contact_robust_cost_coeff*c_covariance;
            fprintf('robust cost function: %4.8f\n',c);
            fprintf('robust cost mean deviation norm-2 value: %4.8f\n',c_mean_dev);
            fprintf('robust cost covariance trace value: %4.8f\n',c_covariance);

            % figure(1)
            % kkk = 2;
            % plot(x(kkk,:),'b-');
            % hold on;
            % for j = 1:n_sampling_point
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
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px(:,:,1) = obj.options.Px_coeff*eye(nx);
                Px_init = Px(:,:,1);
                w_noise = [];
                Pw = [];
                
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
                c = trace(Px_init);
                %c = 1/2*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(nx))) + nx/2*log(2*pi);
                c_covariance = 0;
                c_mean_dev = 0;
                dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(timestep,Sig,u_fdb_k)
                    [xdot,dxdot] = obj.plant.dynamics(0,Sig,u_fdb_k);
                    xdn = Sig + timestep*xdot;
                    df = [xdot,eye(4),zeros(4,nu)]+timestep*full(dxdot);
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
                            c = c + kappa*trace(Px_init);
                            c_covariance = c_covariance + kappa*trace(Px_init);
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
                        %                     if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                        %                         obj.plant.uncertain_mu = w_mu(j);
                        %                     elseif strcmp(obj.plant.uncertainty_source, 'object_initial_position')
                        %                         obj.plant.uncertain_phi = w_phi(:,j);
                        %                         Sig_init(9:10,j,k) = Sig_init(9:10,j,k) + obj.plant.uncertain_phi;%object x and y position uncertainty
                        %                     elseif strcmp(obj.plant.uncertainty_source, 'object_initial_orientation')
                        %                         obj.plant.uncertain_ori = w_ori(:,j);
                        %                         Sig_init(12:13,j,k) = Sig_init(12:13,j,k) + obj.plant.uncertain_ori;%object x and y orientation uncertainty
                        %                     elseif strcmp(obj.plant.uncertainty_source, 'friction_coeff+object_initial_position')
                        %                         obj.plant.uncertain_mu = w_mu(j);
                        %                         obj.plant.uncertain_phi = w_phi(:,j);
                        %                         Sig_init(9:10,j,k) = Sig_init(9:10,j,k) + obj.plant.uncertain_phi;%object x and y position uncertainty
                        %                     elseif strcmp(obj.plant.uncertainty_source, 'friction_coeff+object_initial_orientation')
                        %                         obj.plant.uncertain_mu = w_mu(j);
                        %                         obj.plant.uncertain_ori = w_ori(:,j);
                        %                         Sig_init(12:13,j,k) = Sig_init(12:13,j,k) + obj.plant.uncertain_ori;%object x and y orientation uncertainty
                        %                     end
                        
                        %[H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:nx/2,j,k),Sig_init(nx/2+1:nx,j,k));
                        %Hinv(:,:,j,k) = inv(H);
                        
                        % add feedback control
                        u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                        
                        external_force_index = j;
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
                        %df(:,:,j) = df_numeric;
                        
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
                    c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))) ...
                    %    + nx/2*log(2*pi);%ML mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(nx))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                    %c = c + (x(:,k+1)-x_mean(:,k+1))'*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                    c = c + kappa*trace(Px(:,:,k+1));%ML mean deviation version
                    % note that: 1/2 coefficient is removed.
                    c_mean_dev = c_mean_dev + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;
                    c_covariance = c_covariance + kappa*trace(Px(:,:,k+1));
                    
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
                            dSig_i_kplus1_dx = zeros(nx);
                            dSig_i_kplus1_du = zeros(nx,nu);
                            chain_rule_indx = k-j;
                            if j ~= 1
                                %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                            else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2);%[proved]
                            end
                            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            
                            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
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
                        dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                        dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                        %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                        %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
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
            
        function obj = addRunningCost(obj,running_cost_function)
            % Adds an integrated cost to all time steps, which is
            % numerical implementation specific (thus abstract)
            % this cost is assumed to be time-invariant
            % @param running_cost_function a function handle
            %  of the form running_cost_function(dt,x,u)
            
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            
            for i=1:obj.N-1,
                switch obj.options.integration_method
                    case RobustDirtranTrajectoryOptimization.FORWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)};
                    case RobustDirtranTrajectoryOptimization.BACKWARD_EULER
                        running_cost = FunctionHandleObjective(1+nX+nU, running_cost_function);
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i+1);obj.u_inds(:,i)};
                    case RobustDirtranTrajectoryOptimization.MIDPOINT
                        running_cost = FunctionHandleObjective(1+2*nX+2*nU,...
                            @(h,x0,x1,u0,u1) obj.midpoint_running_fun(running_cost_function,h,x0,x1,u0,u1));
                        inds_i = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.u_inds(:,i+1)};
                    otherwise
                        error('Drake:RobustDirtranTrajectoryOptimization:InvalidArgument','Unknown integration method');
                end
                
                obj = obj.addCost(running_cost,inds_i);
            end
        end
        
        function obj = addRobustRunningCost(obj)
            nX = obj.nX;
            nU = obj.nU;
            N = obj.N;
            
            dim = N-1 + N*nX;%N-1 + N*nX + (N-1)*nU;
            cost = FunctionHandleObjective(dim,@obj.robust_cost,1);
            cost.grad_method = 'user';
            obj = obj.addCost(cost, {obj.h_inds'; obj.x_inds}); % to be double checked
        end
        
        function [c, dc] = robust_cost(obj,r,r_average)
            nX = obj.nX;
            nU = obj.nU;
            nW = obj.nW;
            N = obj.N;
            
            c = 0;
            dc = zeros(1,(N-1)*(1+nX+nU)+nX);
            for k = 1:(N-1)
                c = c + 0.5* (r(:,i) - r_average(:,i))^2;
                
                dc = dc + dcdE*dE(:,:,k) + dcdK*dK(:,:,k);
            end
            c = c + trace(obj.Qrf*E(:,:,N));
            dcdE = vec(obj.Qrf)';
            dc = dc + dcdE*dE(:,:,N);
        end
        
    end
    
    
    methods (Access=protected)
        function [f,df] = forward_constraint_fun(obj,h,x0,x1,u)
            nX = obj.plant.getNumStates();
            [xdot,dxdot] = obj.plant.dynamics(0,x0,u);
            f = x1 - x0 - h*xdot;
            df = [-xdot (-eye(nX) - h*dxdot(:,2:1+nX)) eye(nX) -h*dxdot(:,nX+2:end)];
        end
        
        function [f,df] = backward_constraint_fun(obj,h,x0,x1,u)
            nX = obj.plant.getNumStates();
            [xdot,dxdot] = obj.plant.dynamics(0,x1,u);
            f = x1 - x0 - h*xdot;
            df = [-xdot -eye(nX) (eye(nX) - h*dxdot(:,2:1+nX)) -h*dxdot(:,nX+2:end)];
        end
        
        function [f,df] = midpoint_constraint_fun(obj,h,x0,x1,u0,u1)
            nX = obj.plant.getNumStates();
            [xdot,dxdot] = obj.plant.dynamics(0,.5*(x0+x1),.5*(u0+u1));
            f = x1 - x0 - h*xdot;
            df = [-xdot (-eye(nX) - .5*h*dxdot(:,2:1+nX)) (eye(nX)- .5*h*dxdot(:,2:1+nX)) -.5*h*dxdot(:,nX+2:end) -.5*h*dxdot(:,nX+2:end)];
        end
        
        function [f,df] = midpoint_running_fun(obj,running_handle,h,x0,x1,u0,u1)
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            [f,dg] = running_handle(h,.5*(x0+x1),.5*(u0+u1));
            
            df = [dg(:,1) .5*dg(:,2:1+nX) .5*dg(:,2:1+nX) .5*dg(:,2+nX:1+nX+nU) .5*dg(:,2+nX:1+nX+nU)];
        end
    end
end