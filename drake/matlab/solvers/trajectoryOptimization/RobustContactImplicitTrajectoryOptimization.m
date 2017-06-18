classdef RobustContactImplicitTrajectoryOptimization < DirectTrajectoryOptimization
    % phi, lambda
    properties
        nC
        nD % number of friction elements per contact
        
        l_inds % orderered [lambda_N;lambda_f1;lambda_f2;...;gamma] for each contact sequentially
        lfi_inds % nD x nC indexes into lambda for each time step
        LCP_slack_inds % slack variable for LCP component <z,f(x,z) = LCP_slack > = 0
        lambda_mult
        ljl_inds  % joint limit forces
        jl_lb_ind  % joint indices where the lower bound is finite
        jl_ub_ind % joint indices where the lower bound is finite
        nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
        
        nonlincompl_constraints
        nonlincompl_slack_inds
        
        nonlincompl_constraints_purturb
        nonlincompl_slack_inds_purturb
        
        % complementarity matrix
        W
        r
        M
    end
    
    properties (Constant)
        FORWARD_EULER = 1;
        BACKWARD_EULER = 2;
        MIDPOINT = 3;  % DEFAULT
        MIXED = 4;   % matched to TimeSteppingRigidBodyManipulator. Forward on qd, backward on q
    end
    
    methods
        function obj = RobustContactImplicitTrajectoryOptimization(plant,N,duration,options)
            if nargin<4, options=struct(); end
            
            if ~isfield(options,'nlcc_mode')
                options.nlcc_mode = 2;
            end
            if ~isfield(options,'lincc_mode')
                options.lincc_mode = 1;
            end
            if ~isfield(options,'compl_slack')
                options.compl_slack = 0;
            end
            if ~isfield(options,'lincompl_slack')
                options.lincompl_slack = 0;
            end
            if ~isfield(options,'jlcompl_slack')
                options.jlcompl_slack = 0;
            end
            if ~isfield(options,'lambda_mult')
                options.lambda_mult = 1;
            end
            if ~isfield(options,'lambda_jl_mult')
                options.lambda_jl_mult = 1;
            end
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = RobustContactImplicitTrajectoryOptimization.MIDPOINT;
            end
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
            
        end
        
        function obj = addDynamicConstraints(obj)
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            nq = obj.plant.getNumPositions();
            N = obj.N;
            
            constraints = cell(N-1,1);
            foot_horizontal_distance_constraints = cell(N-1,1);
            foot_height_diff_constraints = cell(N-1,1);
            lincompl_constraints = cell(N-1,1);
            obj.nonlincompl_constraints = cell(N-1,1);
            obj.nonlincompl_slack_inds = cell(N-1,1);
            jlcompl_constraints = cell(N-1,1);
            dyn_inds = cell(N-1,1);
            
            n_vars = 2*nX + nU + 1 + obj.nC*(2+obj.nD) + obj.nJL;
            cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.dynamics_constraint_fun);
            q0 = getZeroConfiguration(obj.plant);
            
            cnstr_foot_horizontal_distance = FunctionHandleConstraint(-inf(2,1),zeros(2,1),nX,@obj.foot_horizontal_distance_constraint_fun);
            cnstr_foot_height_diff = FunctionHandleConstraint(-inf(2,1),zeros(2,1),nX,@obj.foot_height_diff_constraint_fun);
            cnstr_CoM_vertical_velocity = FunctionHandleConstraint(-inf(1,1),zeros(1,1),1,@obj.CoM_vertical_velocity_fun);
            
            [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(q0,false,obj.options.active_collision_options);
            
            for i=1:obj.N-1,
                dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.l_inds(:,i);obj.ljl_inds(:,i)};
                constraints{i} = cnstr;
                obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                
                % add foot horizontal distance constraint
                foot_horizontal_distance_inds{i} = {obj.x_inds(:,i)};
                foot_horizontal_distance_constraints{i} = cnstr_foot_horizontal_distance;
                obj = obj.addConstraint(foot_horizontal_distance_constraints{i}, foot_horizontal_distance_inds{i});
                
                % add foot height diff constraint
                foot_height_diff_inds{i} = {obj.x_inds(:,i)};
                foot_height_diff_constraints{i} = cnstr_foot_height_diff;
                obj = obj.addConstraint(foot_height_diff_constraints{i}, foot_height_diff_inds{i});
                
                % add max CoM vertical velocity constraint
                CoM_vertical_velocity_inds{i} = {obj.x_inds(8,i)};
                CoM_vertical_velocity_constraints{i} = cnstr_CoM_vertical_velocity;
                obj = obj.addConstraint(CoM_vertical_velocity_constraints{i}, CoM_vertical_velocity_inds{i});
                
                if obj.nC > 0
                    % indices for (i) gamma
                    gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end,i);
                    % awkward way to pull out these indices, for (i) lambda_N and
                    % lambda_f
                    lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),i);
                    
                    obj.options.nlcc_mode = 5;% robust mode
                    obj.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode);
                    obj.nonlincompl_slack_inds{i} = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraints{i}.n_slack; % index the six slack variables: gamma in NonlinearComplementarityConstraint
                    obj = obj.addConstraint(obj.nonlincompl_constraints{i},[obj.x_inds(:,i+1);gamma_inds;lambda_inds;obj.LCP_slack_inds(:,i)]);
                    
                    % add ERM cost for sliding velocity constraint uncertainty
                    obj = obj.addCost(FunctionHandleObjective(6+2,@(lambda,gamma)ERMcost_velocity(obj,lambda,gamma),1),{lambda_inds;gamma_inds});
                    
                    % linear complementarity constraint
                    %   gamma /perp mu*lambda_N - sum(lambda_fi)
                    %
                    %  Generate terms W,r,M,gamma_inds so that
                    %  gamma = y(gamma_inds)
                    %  Wz+Mx+r = mu*lambda_N - sum(lambda_fi)
                    r = zeros(obj.nC,1);
                    W = zeros(obj.nC,obj.nC);
                    M = zeros(obj.nC,obj.nC*(1+obj.nD));
                    for k=1:obj.nC,
                        M(k,1 + (k-1)*(1+obj.nD)) = mu(k);
                        M(k,(2:obj.nD+1) + (k-1)*(1+obj.nD)) = -ones(obj.nD,1);
                    end
                    
                    %lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode);
                    %obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds;obj.LCP_slack_inds(:,i)]);
                    
                    % add expected residual minimization cost
                    obj.W = W;
                    obj.r = r;
                    obj.M = M;
                    
                    assert(size(lambda_inds,1) == 6);
                    assert(size(gamma_inds,1) == 2);
                    
                    % add ERM cost for friction cone coefficient uncertainty
                    obj = obj.addCost(FunctionHandleObjective(6+2,@(lambda,gamma)ERMcost_friction(obj,lambda,gamma),1),{lambda_inds;gamma_inds});
                    
                end
                
                if obj.nJL > 0
                    % joint limit linear complementarity constraint
                    % lambda_jl /perp [q - lb_jl; -q + ub_jl]
                    W_jl = zeros(obj.nJL);
                    [r_jl,M_jl] = jointLimitConstraints(obj.plant,q0);
                    jlcompl_constraints{i} = LinearComplementarityConstraint(W_jl,r_jl,M_jl,obj.options.lincc_mode);
                    
                    obj = obj.addConstraint(jlcompl_constraints{i},[obj.x_inds(1:nq,i+1);obj.ljl_inds(:,i);obj.LCP_slack_inds(:,i)]);
                end
            end
            
            if (obj.nC > 0)
                obj = obj.addCost(FunctionHandleObjective(length(obj.LCP_slack_inds),@(slack)robustLCPcost(obj,slack),1),obj.LCP_slack_inds(:));
            end
            
            function [f,df] = ERMcost_friction(obj, lambda, gamma)
                %                 x = y(1:xdim);
                %                 z = y(xdim+1:end-1);
                %                 slack_var = y(end);
                
                zdim = size(obj.W,2);
                xdim = size(obj.M,2);
                
                g = obj.W*gamma + obj.M*lambda + obj.r;
                dg = [obj.M obj.W];
                
                % distribution-free version
                delta = 100;% coefficient
                f = delta/2 * norm(gamma.*g)^2;% - slack_var*ones(zdim,1);
                df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
                % probabilistic version
                sigma = 0.5;
                mu_cov = sigma^2*[lambda(1);lambda(4)];
                delta = 100;% coefficient
                f = delta/2 * norm(gamma.*g - mu_cov)^2;% - slack_var*ones(zdim,1);
                df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
            end
            
            function [f,df] = ERMcost_velocity(obj, lambda, gamma)
                %                 x = y(1:xdim);
                %                 z = y(xdim+1:end-1);
                %                 slack_var = y(end);
                
                zdim = size(obj.W,2);
                xdim = size(obj.M,2);
                
                g = obj.W*gamma + obj.M*lambda + obj.r;
                dg = [obj.M obj.W];
                
                % von Mises-Fisher distribution for quaternion rotation vector
                mu_dirc = [1,0,0,0]'; % this is for flat terrain, no rotation. can be terrain dependent.
                kappa = 10;
                I_kappa_plus = exp(kappa) + exp(-kappa);
                I_kappa_minus = exp(kappa) - exp(-kappa);
                
                h_kappa = (kappa*I_kappa_plus - I_kappa_minus)/(kappa*I_kappa_minus);
                mu_r = mu_dirc*h_kappa;
                Sigma_r = h_kappa/kappa * eye(4) + (1 - 3*h_kappa/kappa - h_kappa^2)*mu_q*mu_q';
                
                mu_w = mu_r(1);mu_x = mu_r(2);mu_y = mu_r(3);mu_z = mu_r(4);
                
                Rbar = [1-2*mu_y^2-2*mu_z^2, 2*mu_x*mu_y-2*mu_z*mu_w, 2*mu_x*mu_z+2*mu_y*mu_w;
                    2*mu_x*mu_y+2*mu_z*mu_w, 1-2*mu_x^2-2*mu_z^2, 2*mu_y*mu_z-2*mu_x*mu_w;
                    2*mu_x*mu_z-2*mu_y*mu_w, 2*mu_y*mu_z+2*mu_x*mu_w, 1-2*mu_x^2-2*mu_y^2];
                
                % normal direction n
                F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
                    -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
                    2*mu_w, -2*mu_x, -2*mu_y, 2*mu_z];
                
                Fc = Rbar(:,3) - F*mu_q;
                
                % tangential direction D_{r,x}
                G = [2*mu_w, 2*mu_x, -2*mu_y, -2*mu_z;
                    2*mu_z, 2*mu_y, 2*mu_x,  2*mu_w;
                    -2*mu_y, 2*mu_z, -2*mu_w, 2*mu_x];
                
                Gc = Rbar(:,1) - G*mu_q;
                
                % tangential direction D_{r,y}
                H = [-2*mu_z, 2*mu_y,  2*mu_x, -2*mu_w;
                    2*mu_w,  -2*mu_x, 2*mu_y, -2*mu_z;
                    2*mu_x,  2*mu_w,  2*mu_z, 2*mu_y];
                
                Hc = Rbar(:,2) - H*mu_q;
                
                % Manipulator dynamics
                % ...
                Minv = eye(6);%inv(M);
                for i=1:6
                    dM(:,:,i) = eye(6);
                end
                Jg = ones(3,6);
                
                
                % Composed matrices
                K = Minv*Jg'*G;
                L = Minv*Jg'*Gc;
                Wx = K'*Jg'*G;
                Yx = K'*Jg'*Gc;
                
                % expectation of M_v_x
                
                E_Drx_Drx = trace(Wx*Sigma_r) + mu_r'*Wx*mu_r + L'*Jg'*(2*G*mu_r + Gc);
                V_Drx_Drx = 2*trace(K*Sigma_r*Wx*Sigma_r*G'*Jg) + 4*(mu_r'*Wx + Yx')*Sigma_r*(Wx*mu_r + Yx) ...
                    +(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc)^2 - E_Drx_Drx^2;
                
                for i=1:6
                    dE_Drx_Drx(i) = trace( (-Minv*Jg'*(G*(Sigma_r + mu_r'*mu_r)*G' + Gc*(2*G*mu_r + Gc)')*Jg*Minv)'*dM(:,:,i) );
                    dV_Drx_Drx_first_chain(:,:,i) = -4*K*Sigma_r*Wx*Sigma_r*K' + 4*(-K*mu_r*mu_r'*Wx*Sigma_r*K' - K*Sigma_r*Wx*mu_r*mu_r'*K' ...
                        -2*K*mu_r*Yx'*Sigma_r*K'-2*K*Sigma_r*Wx*mu_r*L'-L*Yx'*Sigma_r*K'-K*Sigma_r*Yx*L') ...
                        +2*(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc) ...
                        *(-K*(Sigma_r + mu_r*mu_r')*K' - 2*K*mu_r*L' - L*L');
                    dV_Drx_Drx(i) = trace(dV_Drx_Drx_first_chain(:,:,i)'*dM(:,:,i));
                    dV_Drx_Drx(i) = dV_Drx_Drx(i) - 2*E_Drx_Drx*dE_Drx_Drx(i);
                end
                
                % be careful with constant h
                
                % distribution-free version
                delta = 100;% coefficient
                f = delta/2 * norm(gamma.*g)^2;% - slack_var*ones(zdim,1);
                df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
                % probabilistic version
                sigma = 0.5;
                mu_cov = sigma^2*[lambda(1);lambda(4)];
                delta = 100;% coefficient
                f = delta/2 * norm(gamma.*g - mu_cov)^2;% - slack_var*ones(zdim,1);
                df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
            end
            
            % nonlinear complementarity constraints:
            %   lambda_N /perp phi(q)
            %   lambda_fi /perp gamma + Di*psi(q,v)
            % x = [q;v;gamma]
            % z = [lambda_N;lambda_F1;lambda_f2] (each contact sequentially)
            function [f,df] = nonlincompl_fun(y)
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                x = y(1:nq+nv+obj.nC);
                z = y(nq+nv+obj.nC+1:end);
                gamma = x(nq+nv+1:end);
                q = x(1:nq);
                v = x(nq+1:nq+nv);
                
                [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                
                f = zeros(obj.nC*(1+obj.nD),1);
                df = zeros(obj.nC*(1+obj.nD),nq+nv+obj.nC*(2+obj.nD));
                
                f(1:1+obj.nD:end) = phi;
                df(1:1+obj.nD:end,1:nq) = n;
                for j=1:obj.nD,
                    f(1+j:1+obj.nD:end) = gamma+D{j}*v;
                    df(1+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                    df(1+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                    df(1+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
                end
            end
        end
        
        function [c,dc] = robustLCPcost(obj, slack_var)
            c = 1000*sum(slack_var);
            dc = 1000*ones(1,length(slack_var));
        end
        
        function [f,df] = dynamics_constraint_fun(obj,h,x0,x1,u,lambda,lambda_jl)
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nu = obj.plant.getNumInputs;
            nl = length(lambda);
            njl = length(lambda_jl);
            
            lambda = lambda*obj.options.lambda_mult;
            lambda_jl = lambda_jl*obj.options.lambda_jl_mult;
            
            assert(nq == nv) % not quite ready for the alternative
            
            q0 = x0(1:nq);
            v0 = x0(nq+1:nq+nv);
            q1 = x1(1:nq);
            v1 = x1(nq+1:nq+nv);
            
            switch obj.options.integration_method
                case RobustContactImplicitTrajectoryOptimization.MIDPOINT
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics((q0+q1)/2,(v0+v1)/2);
                    dH0 = dH/2;
                    dC0 = dC/2;
                    dB0 = dB/2;
                    dH1 = dH/2;
                    dC1 = dC/2;
                    dB1 = dB/2;
                case RobustContactImplicitTrajectoryOptimization.FORWARD_EULER
                    [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                    dH1 = zeros(nq^2,2*nq);
                    dC1 = zeros(nq,2*nq);
                    dB1 = zeros(nq*nu,2*nq);
                case RobustContactImplicitTrajectoryOptimization.BACKWARD_EULER
                    [H,C,B,dH1,dC1,dB1] = obj.plant.manipulatorDynamics(q1,v1);
                    dH0 = zeros(nq^2,2*nq);
                    dC0 = zeros(nq,2*nq);
                    dB0 = zeros(nq*nu,2*nq);
                case RobustContactImplicitTrajectoryOptimization.MIXED
                    [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                    dH1 = zeros(nq^2,2*nq);
                    dC1 = zeros(nq,2*nq);
                    dB1 = zeros(nq*nu,2*nq);
            end
            
            BuminusC = B*u-C;
            if nu>0,
                dBuminusC0 = matGradMult(dB0,u) - dC0;
                dBuminusC1 = matGradMult(dB1,u) - dC1;
            else
                dBuminusC0 = -dC0;
                dBuminusC1 = -dC1;
            end
            
            switch obj.options.integration_method
                case RobustContactImplicitTrajectoryOptimization.MIDPOINT
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*(v0 + v1)/2;
                    dfq = [-(v1+v0)/2, -eye(nq), -h/2*eye(nq), eye(nq), -h/2*eye(nq) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization.FORWARD_EULER
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*v0;
                    dfq = [-v0, -eye(nq), -h*eye(nq), eye(nq), zeros(nq,nv) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization.BACKWARD_EULER
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*v1;
                    dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization.MIXED
                    fq = q1 - q0 - h*v1;
                    dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq) zeros(nq,nu+nl+njl)];
            end
            
            
            % H*v1 = H*v0 + h*(B*u - C) + n^T lambda_N + d^T * lambda_f
            fv = H*(v1 - v0) - h*BuminusC;
            % [h q0 v0 q1 v1 u l ljl]
            
            dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl)] + ...
                [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,nu+nl+njl)];
            
            if nl>0
                [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                % construct J and dJ from n,D,dn, and dD so they relate to the
                % lambda vector
                J = zeros(nl,nq);
                J(1:2+obj.nD:end,:) = n;
                dJ = zeros(nl*nq,nq);
                dJ(1:2+obj.nD:end,:) = dn;
                
                for j=1:length(D),
                    J(1+j:2+obj.nD:end,:) = D{j};
                    dJ(1+j:2+obj.nD:end,:) = dD{j};
                end
                
                fv = fv - J'*lambda;
                dfv(:,2+nq+nv:1+2*nq+nv) = dfv(:,2+nq+nv:1+2*nq+nv) - matGradMult(dJ,lambda,true);
                dfv(:,2+2*nq+2*nv+nu:1+2*nq+2*nv+nu+nl) = -J'*obj.options.lambda_mult;
            end
            
            if njl>0
                [~,J_jl] = jointLimitConstraints(obj.plant,q1);
                
                fv = fv - J_jl'*lambda_jl;
                dfv(:,2+2*nq+2*nv+nu+nl:1+2*nq+2*nv+nu+nl+njl) = -J_jl'*obj.options.lambda_jl_mult;
            end
            
            f = [fq;fv];
            df = [dfq;dfv];
        end
        
        function [f,df] = foot_horizontal_distance_constraint_fun(obj,x)
            nv = obj.plant.getNumVelocities;
            
            % hard coding fwd kinematics
            CoM_x_pos = x(1);
            q_stance_hip = x(3);
            q_stance_knee = x(4);
            q_swing_hip = x(5);
            q_swing_knee = x(6);
            l_thigh  = 0.5;
            l_calf  = 0.5;
            
            %swing foot sagittal position
            x_swing = CoM_x_pos - l_thigh*sin(q_swing_hip) - l_calf*sin(q_swing_knee + q_swing_hip);
            x_stance = CoM_x_pos - l_thigh*sin(q_stance_hip) - l_calf*sin(q_stance_hip + q_stance_knee);
            
            foot_horizontal_distance_max = 0.2;
            CoM_swing_foot_horizontal_distance_max = 0.4;
            CoM_stance_foot_horizontal_distance_max = 0.2;
            
            f = [x_swing - x_stance - foot_horizontal_distance_max;
                x_swing - CoM_x_pos - CoM_swing_foot_horizontal_distance_max];
            % x_stance - CoM_x_pos - CoM_stance_foot_horizontal_distance_max
            
            df = [zeros(1,2), l_thigh*cos(q_stance_hip) + l_calf*cos(q_stance_hip + q_stance_knee), l_calf*cos(q_stance_hip + q_stance_knee), ...
                -l_thigh*cos(q_swing_hip) - l_calf*cos(q_swing_knee + q_swing_hip), -l_calf*cos(q_stance_hip + q_stance_knee), zeros(1,nv);
                zeros(1,4), -l_thigh*cos(q_swing_hip) - l_calf*cos(q_swing_knee + q_swing_hip), -l_calf*cos(q_stance_hip + q_stance_knee), zeros(1,nv)];
            % zeros(1,2), -l_thigh*cos(q_stance_hip) - l_calf*cos(q_stance_hip + q_stance_knee), -l_calf*cos(q_stance_hip + q_stance_knee), zeros(1,2+nv)
        end
        
        function [f,df] = foot_height_diff_constraint_fun(obj,x)
            nv = obj.plant.getNumVelocities;
            
            % hard coding fwd kinematics
            CoM_z_pos = x(1);
            q_stance_hip = x(3);
            q_stance_knee = x(4);
            q_swing_hip = x(5);
            q_swing_knee = x(6);
            l_thigh  = 0.5;
            l_calf  = 0.5;
            
            %swing foot vertical position
            z_swing = CoM_z_pos - l_thigh*cos(q_swing_hip) - l_calf*cos(q_swing_knee + q_swing_hip);
            z_stance = CoM_z_pos - l_thigh*cos(q_stance_hip) - l_calf*cos(q_stance_hip + q_stance_knee);
            
            foot_height_distance_max = 0.5;
            CoM_foot_height_diff_max = 0.6;
            
            f = [z_swing - z_stance - foot_height_distance_max;
                CoM_foot_height_diff_max - CoM_z_pos + z_swing];
            
            df = [zeros(1,2), -l_thigh*sin(q_stance_hip) - l_calf*sin(q_stance_hip + q_stance_knee), -l_calf*sin(q_stance_hip + q_stance_knee), ...
                l_thigh*sin(q_swing_hip) + l_calf*sin(q_swing_knee + q_swing_hip), l_calf*sin(q_stance_hip + q_stance_knee), zeros(1,nv);
                zeros(1,4), l_thigh*sin(q_swing_hip) + l_calf*sin(q_swing_knee + q_swing_hip), l_calf*sin(q_stance_hip + q_stance_knee), zeros(1,nv)];
        end
        
        function [f,df] = CoM_vertical_velocity_fun(obj,CoM_z_vel)
            CoM_vertical_velocity_max = 0.4;
            
            f = [CoM_z_vel];
            df = [ones(1)];
        end
        
        function [xtraj,utraj,ltraj,ljltraj,slacktraj,z,F,info,infeasible_constraint_name] = solveTraj(obj,t_init,traj_init)
            [xtraj,utraj,z,F,info,infeasible_constraint_name] = solveTraj@DirectTrajectoryOptimization(obj,t_init,traj_init);
            t = [0; cumsum(z(obj.h_inds))];
            if obj.nC>0
                ltraj = PPTrajectory(foh(t,reshape(z(obj.l_inds),[],obj.N)));
                slacktraj = PPTrajectory(foh(t,[reshape(z(obj.LCP_slack_inds),size(obj.LCP_slack_inds,1),obj.N-1),z(obj.LCP_slack_inds(:,end))]));
            else
                ltraj = [];
                slacktraj = [];
            end
            if obj.nJL>0
                ljltraj = PPTrajectory(foh(t,reshape(z(obj.ljl_inds),[],obj.N)));
            else
                ljltraj = [];
            end
        end
        
        function obj = setupVariables(obj,N)
            obj = setupVariables@DirectTrajectoryOptimization(obj,N);
            [~,normal,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), false, obj.options.active_collision_options);
            obj.nC = size(normal, 2);
            obj.nD = 2*length(d);
            
            nContactForces = obj.nC*(2 + obj.nD);
            
            obj.l_inds = reshape(obj.num_vars + (1:N * nContactForces),nContactForces,N);
            obj = obj.addDecisionVariable(N * nContactForces);
            
            obj.lfi_inds = zeros(obj.nD,obj.nC);
            for i=1:obj.nC,
                obj.lfi_inds(:,i) = (2:1+obj.nD)' + (i-1)*(2+obj.nD)*ones(obj.nD,1);
            end
            
            obj.nJL = obj.plant.getNumJointLimitConstraints();
            obj.ljl_inds = reshape(obj.num_vars + (1:N * obj.nJL),obj.nJL,N);
            
            % joint limit constraints
            [jl_lb,jl_ub] = obj.plant.getJointLimits();
            obj.jl_lb_ind = find(jl_lb ~= -inf);
            obj.jl_ub_ind = find(jl_ub ~= inf);
            
            obj = obj.addDecisionVariable(N * obj.nJL);
            
            %add maximal LCP slack variable into decision variable
            %obj.LCP_slack_inds = reshape(obj.num_vars + (1:2*(N-1)), 2, N-1);%[diff slack var]
            obj.LCP_slack_inds = reshape(obj.num_vars + (1:N-1), 1, N-1);%[single slack var]
            %obj = obj.addDecisionVariable(2*(N-1));%[diff slack var]
            obj = obj.addDecisionVariable(N-1);%[single slack var]
        end
        
        % evaluates the initial trajectories at the sampled times and
        % constructs the nominal z0. Overwrite to implement in a different
        % manner
        function z0 = getInitialVars(obj,t_init,traj_init)
            if isscalar(t_init)
                t_init = linspace(0,t_init,obj.N);
            elseif length(t_init) ~= obj.N
                error('The initial sample times must have the same length as property N')
            end
            z0 = zeros(obj.num_vars,1);
            z0(obj.h_inds) = diff(t_init);
            
            if isfield(traj_init,'u')
                z0(obj.u_inds) = traj_init.u.eval(t_init);
            else
                nU = getNumInputs(obj.plant);
                z0(obj.u_inds) = 0.01*randn(nU,obj.N);
            end
            
            if isfield(traj_init,'x')
                z0(obj.x_inds) = traj_init.x.eval(t_init);
            else
                if ~isfield(traj_init,'u')
                    traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,reshape(z0(obj.u_inds),nU,obj.N))),getInputFrame(obj.plant));
                end
                
                % todo: if x0 and xf are equality constrained, then initialize with
                % a straight line from x0 to xf (this was the previous behavior)
                
                %simulate
                sys_ol = cascade(traj_init.u,obj.plant);
                [~,x_sim] = sys_ol.simulate([t_init(1) t_init(end)]);
                z0(obj.x_inds) = x_sim.eval(t_init);
            end
            
            if obj.nC > 0
                if isfield(traj_init,'l')
                    z0(obj.l_inds) = traj_init.l.eval(t_init);
                else
                    z0(obj.l_inds) = 0;
                end
                
                % initialize LCP slack variables
                if isfield(traj_init,'LCP_slack')
                    LCP_slack_eval = traj_init.LCP_slack.eval(t_init);
                    z0(obj.LCP_slack_inds) = LCP_slack_eval(:,1:obj.N-1); % take N-1 values
                else
                    z0(obj.LCP_slack_inds) = 0;
                end
            end
            if obj.nJL > 0
                if isfield(traj_init,'ljl')
                    z0(obj.ljl_inds) = traj_init.ljl.eval(t_init);
                else
                    z0(obj.ljl_inds) = 0;
                end
            end
            
            if obj.nC > 0
                for i=1:obj.N-1,
                    gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end,i);
                    lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),i);
                    if ~isempty(obj.nonlincompl_slack_inds{i})
                        z0(obj.nonlincompl_slack_inds{i}) = obj.nonlincompl_constraints{i}.slack_fun(z0([obj.x_inds(:,i+1);gamma_inds;lambda_inds]));
                    end
                end
            end
        end
        
        function obj = addRunningCost(obj,running_cost_function)
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            running_cost = FunctionHandleObjective(1+nX+nU,running_cost_function);
            for i=1:obj.N-1,
                obj = obj.addCost(running_cost,{obj.h_inds(i);obj.x_inds(:,i);obj.u_inds(:,i)});
            end
        end
        
        function obj = addRobustLCPConstraints(obj,plant_perturb, perturb_index, SampleNum)
            nX = plant_perturb.getNumStates();
            nU = plant_perturb.getNumInputs();
            nq = plant_perturb.getNumPositions();
            N = obj.N;
            
            lincompl_constraints = cell(N-1,1);
            obj.nonlincompl_constraints_purturb = cell(N-1,1);
            obj.nonlincompl_slack_inds_purturb = cell(N-1,1);
            jlcompl_constraints = cell(N-1,1);
            
            q0 = getZeroConfiguration(plant_perturb);
            [~,~,~,~,~,~,~,mu] = plant_perturb.contactConstraints(q0,false,obj.options.active_collision_options);
            
            for i=1:obj.N-1,
                if obj.nC > 0
                    % indices for (i) gamma
                    gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end,i);
                    % awkward way to pull out these indices, for (i) lambda_N and
                    % lambda_f
                    lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),i);
                    
                    obj.nonlincompl_constraints_purturb{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode,obj.options.compl_slack);
                    obj.nonlincompl_slack_inds_purturb{i} = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraints_purturb{i}.n_slack;
                    obj = obj.addConstraint(obj.nonlincompl_constraints_purturb{i},[obj.x_inds(:,i+1);gamma_inds;lambda_inds]);
                    
                    % linear complementarity constraint
                    %   gamma /perp mu*lambda_N - sum(lambda_fi)
                    %
                    %  Generate terms W,r,M,gamma_inds so that
                    %  gamma = y(gamma_inds)
                    %  Wz+Mx+r = mu*lambda_N - sum(lambda_fi)
                    r = zeros(obj.nC,1);
                    W = zeros(obj.nC,obj.nC);
                    M = zeros(obj.nC,obj.nC*(1+obj.nD));
                    for k=1:obj.nC,
                        M(k,1 + (k-1)*(1+obj.nD)) = mu(k);
                        M(k,(2:obj.nD+1) + (k-1)*(1+obj.nD)) = -ones(obj.nD,1);
                    end
                    
                    lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode,obj.options.lincompl_slack);
                    obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds]);
                end
                
                if obj.nJL > 0
                    % joint limit linear complementarity constraint
                    % lambda_jl /perp [q - lb_jl; -q + ub_jl]
                    W_jl = zeros(obj.nJL);
                    [r_jl,M_jl] = jointLimitConstraints(plant_perturb,q0);
                    jlcompl_constraints{i} = LinearComplementarityConstraint(W_jl,r_jl,M_jl,obj.options.lincc_mode,obj.options.jlcompl_slack);
                    
                    obj = obj.addConstraint(jlcompl_constraints{i},[obj.x_inds(1:nq,i+1);obj.ljl_inds(:,i)]);
                end
            end
            
            % nonlinear complementarity constraints:
            %   lambda_N /perp phi(q)
            %   lambda_fi /perp gamma + Di*psi(q,v)
            % x = [q;v;gamma]
            % z = [lambda_N;lambda_F1;lambda_f2] (each contact sequentially)
            function [f,df] = nonlincompl_fun(y)
                disp('perturbed plant index')
                perturb_index
                
                nq = plant_perturb.getNumPositions;
                nv = plant_perturb.getNumVelocities;
                x = y(1:nq+nv+obj.nC);
                z = y(nq+nv+obj.nC+1:end);
                gamma = x(nq+nv+1:end);
                
                q = x(1:nq);
                v = x(nq+1:nq+nv);
                
                [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = plant_perturb.contactConstraints(q,false,obj.options.active_collision_options);
                
                f = zeros(obj.nC*(1+obj.nD),1);
                df = zeros(obj.nC*(1+obj.nD),nq+nv+obj.nC*(2+obj.nD));
                
                f(1:1+obj.nD:end) = phi;
                df(1:1+obj.nD:end,1:nq) = n;
                for j=1:obj.nD,
                    f(1+j:1+obj.nD:end) = gamma+D{j}*v;
                    df(1+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                    df(1+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                    df(1+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
                end
                f = f./SampleNum;
                df = df./SampleNum;
            end
        end
        
    end
end