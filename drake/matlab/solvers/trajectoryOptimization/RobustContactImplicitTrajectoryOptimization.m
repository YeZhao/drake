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
                    
                    % add expected residual minimization cost
                    obj.W = W;
                    obj.r = r;
                    obj.M = M;
                    
                    % add ERM cost for sliding velocity constraint uncertainty
                    obj = obj.addCost(FunctionHandleObjective(2*nX+nU+6+2+1,@(h,x0,x1,u,lambda,gamma)ERMcost_slidingVelocity(obj,h,x0,x1,u,lambda,gamma),1), ...
                          {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);lambda_inds;gamma_inds});
                    
                    
                    %lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode);
                    %obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds;obj.LCP_slack_inds(:,i)]);
                    
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
                mu_cov = sigma^2*[lambda(1)^2;lambda(4)^2];% make sure it is squared
                delta = 100;% coefficient
                f = delta/2 * (norm(diag(gamma)*g)^2 + norm(mu_cov)^2);% - slack_var*ones(zdim,1);
                df = delta*(diag(gamma)*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]] ...
                    + delta*mu_cov'*[2*sigma^2*lambda(1),zeros(1,7);zeros(1,3), 2*sigma^2*lambda(4), zeros(1,4)];
            end
            
            function [f,df] = ERMcost_slidingVelocity(obj, h, x0, x1, u, lambda, gamma)
                %                 x = y(1:xdim);
                %                 z = y(xdim+1:end-1);
                %                 slack_var = y(end);
                
                zdim = size(obj.W,2);
                xdim = size(obj.M,2);
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nu = obj.plant.getNumInputs;
                nl = length(lambda);
                
                % von Mises-Fisher distribution for quaternion rotation vector
                mu_dirc = [1,0,0,0]'; % this is for flat terrain, no rotation. can be terrain dependent.
                kappa = 10;
                I_kappa_plus = exp(kappa) + exp(-kappa);
                I_kappa_minus = exp(kappa) - exp(-kappa);
                
                h_kappa = (kappa*I_kappa_plus - I_kappa_minus)/(kappa*I_kappa_minus);
                mu_r = mu_dirc*h_kappa;
                Sigma_r = h_kappa/kappa * eye(4) + (1 - 3*h_kappa/kappa - h_kappa^2)*mu_r*mu_r';
                
                mu_w = mu_r(1);mu_x = mu_r(2);mu_y = mu_r(3);mu_z = mu_r(4);
                
                Rbar = [1-2*mu_y^2-2*mu_z^2, 2*mu_x*mu_y-2*mu_z*mu_w, 2*mu_x*mu_z+2*mu_y*mu_w;
                    2*mu_x*mu_y+2*mu_z*mu_w, 1-2*mu_x^2-2*mu_z^2, 2*mu_y*mu_z-2*mu_x*mu_w;
                    2*mu_x*mu_z-2*mu_y*mu_w, 2*mu_y*mu_z+2*mu_x*mu_w, 1-2*mu_x^2-2*mu_y^2];
                
                % normal direction n
                F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
                    -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
                    2*mu_w, -2*mu_x, -2*mu_y, 2*mu_z];
                Fc = Rbar(:,3) - F*mu_r;
                
                % tangential direction D_{r,x}
                G = [2*mu_w, 2*mu_x, -2*mu_y, -2*mu_z;
                    2*mu_z, 2*mu_y, 2*mu_x,  2*mu_w;
                    -2*mu_y, 2*mu_z, -2*mu_w, 2*mu_x];
                Gc = Rbar(:,1) - G*mu_r;
                
                % tangential direction D_{r,y}
                H = [-2*mu_z, 2*mu_y,  2*mu_x, -2*mu_w;
                    2*mu_w,  -2*mu_x, 2*mu_y, -2*mu_z;
                    2*mu_x,  2*mu_w,  2*mu_z, 2*mu_y];
                Hc = Rbar(:,2) - H*mu_r;
                
                % RigidBodyManipulator dynamics
                q0 = x0(1:nq);
                v0 = x0(nq+1:nq+nv);
                q1 = x1(1:nq);
                v1 = x1(nq+1:nq+nv);
                
                % currently only implement MIDPOINT method
                [M,C,B,dM,dC,dB] = obj.plant.manipulatorDynamics((q0+q1)/2,(v0+v1)/2);
                Minv = inv(M);
                
                for i=1:nq
                    dMdq(:,:,i) = reshape(dM(:,i),[nq, nq]);
                    dCdq(:,i) = reshape(dC(:,i),[nq, 1]);
                    dCdqdot(:,i) = reshape(dC(:,i+nq),[nq, 1]);
                end
                
                qdot_prev = v0;
                u_prev = u;
                
                if nl>0
                    [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                    % construct J and dJ from n,D,dn, and dD so they relate to the
                    % lambda vector
                    J = zeros(nl,nq);
                    J(1:1+obj.nD:end,:) = n;
                    dJ = zeros(nl*nq,nq);
                    dJ(1:1+obj.nD:end,:) = dn;
                    
                    for j=1:length(D),
                        J(1+j:1+obj.nD:end,:) = D{j};
                        dJ(1+j:1+obj.nD:end,:) = dD{j};
                    end
                    Jg = J;% update notation
                end
                Jg = J(1:3,:);% temporarily define it for one foot contact point
                
                % Composed matrices
                Fy = H;
                Fyc = Hc;
                U = Minv*Jg'*F;
                V = U'*Jg'*G;
                Vy = U'*Jg'*Fy;
                Vyy = Fy'*Jg*Minv*Jg'*G;
                Z = Minv*Jg'*Fc;
                Zy = Minv*Jg'*Fyc;
                X = G'*Jg*Z;
                Xy = Fy'*Jg*Z;
                Xyy = G'*Jg*Zy;
                K = Minv*Jg'*G;
                Ky = Minv*Jg'*Fy;
                L = Minv*Jg'*Gc;
                Ly = Minv*Jg'*Fyc;
                Wx = K'*Jg'*G;
                Wy = Ky'*Jg'*Fy;
                Yx = K'*Jg'*Gc;
                Yy = Ky'*Jg'*Fyc;
                Q = U'*Jg'*Gc;
                Qy = U'*Jg'*Fyc;
                Qyy = Ky'*Jg'*Gc;
                O = (V+V')*mu_r+X+Q;
                Oy = (Vy+Vy')*mu_r+Xy+Qy;
                Oyy = (Vyy+Vyy')*mu_r+Xyy+Qyy;
                
                J_blk = Jg*(qdot_prev + Minv*(B*u_prev - C)*h);
                T = G'*J_blk;
                Tc = Gc'*J_blk;
                Ty = Fy'*J_blk;
                Tyc = Fyc'*J_blk;
                
                % expectation and covariance of M_v_x
                E_M_Drx_nr = trace(V*Sigma_r) + mu_r'*V*mu_r + Z'*Jg'*(G*mu_r + Gc) + mu_r'*U'*Jg'*Gc;
                V_M_Drx_nr = trace(U*Sigma_r*(V+V')*Sigma_r*G'*Jg) + O'*Sigma_r*O ...
                    +(trace(U*Sigma_r*G'*Jg)+mu_r'*V*mu_r+mu_r'*Q+X'*mu_r+Z'*Jg'*Gc)^2 - E_M_Drx_nr^2;
                
                E_M_Drx_Drx = trace(Wx*Sigma_r) + mu_r'*Wx*mu_r + L'*Jg'*(2*G*mu_r + Gc);
                V_M_Drx_Drx = 2*trace(K*Sigma_r*Wx*Sigma_r*G'*Jg) + 4*(mu_r'*Wx + Yx')*Sigma_r*(Wx*mu_r + Yx) ...
                    +(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc)^2 - E_M_Drx_Drx^2;
                
                E_M_Drx_Dry = trace(Wy*Sigma_r) + mu_r'*Wy*mu_r + Ly'*Jg'*(2*Fy*mu_r + Fyc);
                E_M_Dry_Drx = E_M_Drx_Dry;
                
                for i=1:nq
                    % expectation derivative w.r.t q and qdot
                    dE_M_Drx_nr_dq(i) = trace( (-Minv*Jg'*(F*(Sigma_r + mu_r*mu_r')*G' + Fc*(G*mu_r + Gc)' + F*mu_r*Gc')*Jg*Minv)'*dMdq(:,:,i) );
                    dE_M_Drx_nr_dqdot(i) = 0;
                    % [ToDo: derivative w.r.t Jacobian]
                    
                    dE_M_Drx_Drx_dq(i) = trace( (-Minv*Jg'*(G*(Sigma_r + mu_r*mu_r')*G' + Gc*(2*G*mu_r + Gc)')*Jg*Minv)'*dMdq(:,:,i) );
                    dE_M_Drx_Drx_dqdot(i) = 0;
                    % [ToDo: derivative w.r.t Jacobian]
                    
                    dE_M_Drx_Dry_dq(i) = trace( (-Minv*Jg'*(Fy*(Sigma_r + mu_r*mu_r')*G' + Fyc*(G*mu_r + Gc)' + Fy*mu_r*Gc')*Jg*Minv)'*dMdq(:,:,i) );
                    dE_M_Drx_Dry_dqdot(i) = 0;
                    % [ToDo: derivative w.r.t Jacobian]

                    dE_M_Dry_Drx_dq(i) = dE_M_Drx_Dry_dq(i);
                    dE_M_Dry_Drx_dqdot(i) = 0;
                    
                    % covariance derivative w.r.t q and qdot
                    dV_M_Drx_nr_dq_first_chain(:,:,i) = -K*Sigma_r*(V+V')*Sigma_r*U' - U*Sigma_r*V*Sigma_r*K' - K*Sigma_r*V*Sigma_r*U' ...
                        -U*(mu_r*O'+O*mu_r')*K'-K*(mu_r*O'+O*mu_r')*U'-Z*O'*K'-L*O'*U'-K*O*Z'-U*O*L' ...
                        +2*(trace(U*Sigma_r*G'*Jg)+mu_r'*V*mu_r+mu_r'*Q+X'*mu_r+Z'*Jg'*Gc) ...
                        *(-K*Sigma_r*U' - U*mu_r*mu_r'*K' - U*mu_r*L' - Z*mu_r'*K' - Z*L');
                    dV_M_Drx_nr_dq(i) = trace(dV_M_Drx_nr_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_nr_dq(i) = dV_M_Drx_nr_dq(i) - 2*E_M_Drx_nr*dE_M_Drx_nr_dq(i);%[double check this part]
                    
                    dV_M_Drx_nr_dqdot(i) = 0;
                    
                    dV_M_Drx_Drx_dq_first_chain(:,:,i) = -4*K*Sigma_r*Wx*Sigma_r*K' + 4*(-K*mu_r*mu_r'*Wx*Sigma_r*K' - K*Sigma_r*Wx*mu_r*mu_r'*K' ...
                        -2*K*mu_r*Yx'*Sigma_r*K'-2*K*Sigma_r*Wx*mu_r*L'-L*Yx'*Sigma_r*K'-K*Sigma_r*Yx*L') ...
                        +2*(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc) ...
                        *(-K*(Sigma_r + mu_r*mu_r')*K' - 2*K*mu_r*L' - L*L');
                    dV_M_Drx_Drx_dq(i) = trace(dV_M_Drx_Drx_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_Drx_dq(i) = dV_M_Drx_Drx_dq(i) - 2*E_M_Drx_Drx*dE_M_Drx_Drx_dq(i);%[double check this part]
                    
                    dV_M_Drx_Drx_dqdot(i) = 0;
                    
                    dV_M_Drx_Dry_dq_first_chain(:,:,i) = -K*Sigma_r*(Vyy+Vyy')*Sigma_r*Ky' - Ky*Sigma_r*Vyy*Sigma_r*K' - K*Sigma_r*Vyy*Sigma_r*Ky' ...
                        -Ky*(mu_r*Oyy'+Oyy*mu_r')*K'-K*(mu_r*Oyy'+Oyy*mu_r')*Ky'-Zy*Oyy'*K'-L*Oyy'*Ky'-K*Oyy*Zy'-Ky*Oyy*L' ...
                        +2*(trace(Ky*Sigma_r*G'*Jg)+mu_r'*Vyy*mu_r+mu_r'*Qyy+Xyy'*mu_r+Zy'*Jg'*Gc) ...
                        *(-K*Sigma_r*Ky' - Ky*mu_r*mu_r'*K' - Ky*mu_r*L' - Zy*mu_r'*K' - Zy*L');
                    dV_M_Drx_Dry_dq(i) = trace(dV_M_Drx_Dry_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_Dry_dq(i) = dV_M_Drx_Dry_dq(i) - 2*E_M_Drx_Dry*dE_M_Drx_Dry_dq(i);%[double check this part]
                    
                    dV_M_Drx_Dry_dqdot(i) = 0;
                    
                    dV_M_Dry_Drx_dq(i) = dV_M_Drx_Dry_dq(i);
                    dV_M_Dry_Drx_dqdot(i) = 0;
                end
                
                % expectation and covariance of b_v_x
                E_b_Drx = (mu_r'*G' + Gc')*J_blk;
                V_b_Drx = trace(T*T'*Sigma_r);
                
                for i=1:nq
                    dE_b_Drx_dq(i) = trace( (-h*(K*mu_r + L)*(B*u_prev - C)'*Minv)'*dMdq(:,:,i) ) - trace( h*(K*mu_r + L)'*dCdq(:,i));
                    dE_b_Drx_dqdot(i) = - trace( h*(K*mu_r + L)'*dCdqdot(:,i));
                    
                    dV_b_Drx_dq(i) = trace( (-h*Minv*(B*u_prev - C)*T'*Sigma_r*K' -h*K*Sigma_r*T*(B*u_prev - C)'*Minv)'*dMdq(:,:,i)) ...
                        - trace( (2*h*K*Sigma_r*T)'*dCdq(:,i));
                    dV_b_Drx_dqdot(i) = - trace( (2*h*K*Sigma_r*T)'*dCdqdot(:,i));
                    
                end
                dE_b_Drx_dlambda_n = 0;
                dE_b_Drx_dlambda_t = zeros(2,1);
                dE_b_Drx_dgamma = 0;
                
                dV_b_Drx_dlambda_n = 0;
                dV_b_Drx_dlambda_t = zeros(2,1);
                dV_b_Drx_dgamma = 0;
                
                % be careful with constant h
                
                lambda_n = lambda(1);
                lambda_tx = lambda(2);
                lambda_ty = lambda(3);
                gamma_left = gamma(1);
                lambda_vec = [lambda_n;lambda_tx;lambda_ty;gamma_left];

                E_Phi = zeros(4,1);
                
                E_M_v_x = [h*E_M_Drx_nr, h*E_M_Drx_Drx, h*E_M_Drx_Dry, 1];
                E_Mvx_lambda_plus_bvx = E_M_v_x*lambda_vec + E_b_Drx;
                
                % NCP residual, currently assume no smoothing func applied
                E_Phi(1) = lambda_tx*E_Mvx_lambda_plus_bvx;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dE_Phi_dh(1) = lambda_tx*(E_M_v_x*lambda_vec/h + (mu_r'*G' + Gc')*Jg*Minv*(B*u_prev - C));% the last part is dE_b_Drx/dh
                dE_Phi_dq0(:,1) = lambda_tx*h*(dE_M_Drx_nr_dq*lambda_n+dE_M_Drx_Drx_dq*lambda_tx+dE_M_Drx_Dry_dq*lambda_ty)/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dv0(:,1) = lambda_tx*h*(dE_M_Drx_nr_dqdot*lambda_n+dE_M_Drx_Drx_dqdot*lambda_tx+dE_M_Drx_Dry_dqdot*lambda_ty)/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dq1(:,1) = dE_Phi_dq0(:,1);
                dE_Phi_dv1(:,1) = dE_Phi_dv0(:,1);
                dE_Phi_du(:,1) = lambda_tx*(mu_r'*G' + Gc')*Jg*Minv*B*h;
                dE_Phi_dlambda_n(1) = lambda_tx*h*E_M_Drx_nr;
                dE_Phi_dlambda_tx(1) = E_Mvx_lambda_plus_bvx + h*E_M_Drx_Drx*lambda_tx;
                dE_Phi_dlambda_ty(1) = lambda_tx*h*E_M_Drx_Dry;
                
                dE_Phi(:,1) = [dE_Phi_dh(1); dE_Phi_dq0(:,1); dE_Phi_dv0(:,1); dE_Phi_dq1(:,1); dE_Phi_dv1(:,1); dE_Phi_du(:,1);dE_Phi_dlambda_n(1); ...
                               dE_Phi_dlambda_tx(1);dE_Phi_dlambda_ty(1);zeros(3,1);lambda_tx;0];
                
                
                %--------------- third LCP condition ---------------%
                % expectation and covariance of M_v_y
                E_M_Dry_nr = trace(Vy*Sigma_r) + mu_r'*Vy*mu_r + Z'*Jg'*(Fy*mu_r + Fyc) + mu_r'*U'*Jg'*Fyc;
                E_M_nr_Dry = trace(Vy'*Sigma_r) + mu_r'*Vy'*mu_r + Ly'*Jg'*(F*mu_r + Fc) + mu_r'*Ky'*Jg'*Fc;%==E_M_Dry_nr
                E_M_Dry_Dry = trace(Wy*Sigma_r) + mu_r'*Wy*mu_r + Ly'*Jg'*(2*Fy*mu_r+Fyc);
                
                for i=1:nq
                    % expectation derivative w.r.t q and qdot
                    dE_M_Dry_nr_dq(i) = trace( (-Minv*Jg'*(F*(Sigma_r + mu_r*mu_r')*Fy' + Fc*(Fy*mu_r + Fyc)' + F*mu_r*Fyc')*Jg*Minv)'*dMdq(:,:,i) );
                    dE_M_Dry_nr_dqdot(i) = 0;
                    % [ToDo: derivative w.r.t Jacobian]
                    
                    dE_M_Dry_Dry_dq(i) = trace( (-Minv*Jg'*(Fy*(Sigma_r + mu_r*mu_r')*Fy' + Fyc*(2*Fy*mu_r + Fyc)')*Jg*Minv)'*dMdq(:,:,i) );
                    dE_M_Dry_Dry_dqdot(i) = 0;
                    % [ToDo: derivative w.r.t Jacobian]
                    
                    % covariance derivative w.r.t q and qdot
                    dV_M_Dry_Dry_dq_first_chain(:,:,i) = -4*Ky*Sigma_r*Wy*Sigma_r*Ky' + 4*(-Ky*mu_r*mu_r'*Wy*Sigma_r*Ky' - Ky*Sigma_r*Wy*mu_r*mu_r'*Ky' ...
                        -2*Ky*mu_r*Yy'*Sigma_r*Ky'-2*Ky*Sigma_r*Wy*mu_r*Ly'-Ly*Yy'*Sigma_r*Ky'-Ky*Sigma_r*Yy*Ly') ...
                        +2*(trace(Ky*Sigma_r*Fy'*Jg)+mu_r'*Wy*mu_r+2*mu_r'*Yy+Ly'*Jg'*Fyc) ...
                        *(-Ky*(Sigma_r + mu_r*mu_r')*Ky' - 2*Ky*mu_r*Ly' - Ly*Ly');
                    dV_M_Dry_Dry_dq(i) = trace(dV_M_Dry_Dry_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Dry_Dry_dq(i) = dV_M_Dry_Dry_dq(i) - 2*E_M_Dry_Dry*dE_M_Dry_Dry_dq(i);%[double check this part]
                    
                    dV_M_Dry_Dry_dqdot(i) = 0;
                    
                    dV_M_Dry_nr_dq_first_chain(:,:,i) = -Ky*Sigma_r*(Vy+Vy')*Sigma_r*U' - U*Sigma_r*Vy*Sigma_r*Ky' - Ky*Sigma_r*Vy*Sigma_r*U' ...
                        -U*(mu_r*Oy'+Oy*mu_r')*Ky'-Ky*(mu_r*Oy'+Oy*mu_r')*U'-Z*Oy'*Ky'-Ly*Oy'*U'-Ky*Oy*Z'-U*Oy*Ly' ...
                        +2*(trace(U*Sigma_r*Fy'*Jg)+mu_r'*Vy*mu_r+mu_r'*Qy+Xy'*mu_r+Z'*Jg'*Fyc) ...
                        *(-Ky*Sigma_r*U' - U*mu_r*mu_r'*Ky' - U*mu_r*Ly' - Z*mu_r'*Ky' - Z*Ly');
                    dV_M_Dry_nr_dq(i) = trace(dV_M_Dry_nr_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Dry_nr_dq(i) = dV_M_Dry_nr_dq(i) - 2*E_M_Dry_nr*dE_M_Dry_nr_dq(i);%[double check this part]
                    
                    dV_M_Dry_nr_dqdot(i) = 0;
                end
                
                % expectation and covariance of b_v_y
                E_b_Dry = (mu_r'*Fy' + Fyc')*J_blk;
                V_b_Dry = trace(Ty*Ty'*Sigma_r);
                
                E_M_v_y = [h*E_M_Dry_nr, h*E_M_Dry_Drx, h*E_M_Dry_Dry, 1];
                E_Mvy_lambda_plus_bvy = E_M_v_y*lambda_vec + E_b_Dry;
                
                % NCP residual, currently assume no smoothing func applied
                E_Phi(2) = lambda_ty*E_Mvy_lambda_plus_bvy;
                %dE_Phi(2) =
                
                % LCP variance matrix of V_Mvx_lambda_plus_bvx
                %fourth order expectation 
                [E_xnxn,dE_xnxn_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xnxx,dE_xnxx_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xnxy,dE_xnxy_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                [E_xxxn,dE_xxxn_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xxxx,dE_xxxx_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xxxy,dE_xxxy_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                [E_xyxn,dE_xyxn_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xyxx,dE_xyxx_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_xyxy,dE_xyxy_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                E_Mvx_lambda_lambda_Mvx = h^2*(E_xnxn*lambda_n^2 + E_xnxx*lambda_n*lambda_tx + E_xnxy*lambda_n*lambda_ty ...
                    + E_xxxn*lambda_tx*lambda_n + E_xxxx*lambda_tx^2 + E_xxxy*lambda_tx*lambda_ty ...
                    + E_xyxn*lambda_ty*lambda_n + E_xyxx*lambda_ty*lambda_tx + E_xyxy*lambda_ty^2);
                E_Mvx_lambda_lambda_Mvx = E_Mvx_lambda_lambda_Mvx + 2*h*E_M_Drx_nr*lambda_n*gamma_left + 2*h*E_M_Drx_Drx*lambda_tx*gamma_left ...
                                          + 2*h*E_M_Drx_Dry*lambda_ty*gamma_left + gamma_left^2;
                
                %cov(M_{v,x}^T,b_{v,x}) = E[(M_{v,x}^T - E(M_{v,x}^T))(b_{v,x}-E(b_{v,x}))] = E[M_{v,x}^T*b_{v,x}] - E[M_{v,x}^T]*E[b_{v,x}];
                % E[M_{v,x}^T*b_{v,x}]
                [E_Mvx_bvx_Drx_nr_Drx,dE_Mvx_bvx_Drx_nr_Drx_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(G,Gc,F,Fc,G,Gc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                [E_Mvx_bvx_Drx_Drx_Drx,dE_Mvx_bvx_Drx_Drx_Drx_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(G,Gc,G,Gc,G,Gc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                [E_Mvx_bvx_Drx_Dry_Drx,dE_Mvx_bvx_Drx_Dry_Drx_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(G,Gc,Fy,Fyc,G,Gc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                E_Mvx_bvx = [E_Mvx_bvx_Drx_nr_Drx;E_Mvx_bvx_Drx_Drx_Drx;E_Mvx_bvx_Drx_Dry_Drx;E_b_Drx];
                
                %  V[M_{v,x}*lambda+b_{v,x}] = V[M_{v,x}*lambda] + cov(M_{v,x}*lambda,b_{v,x}) + cov(b_{v,x},M_{v,x}*lambda) + V[b_{v,x}]
                % = V[M_{v,x}*lambda] + 2*lambda^T*cov(M_{v,x}^T,b_{v,x}) + V[b_{v,x}]
                % =E[M_{v,x}*lambda*lambda^T*M_{v,x}^T] - (E[M_{v,x}*lambda])^2 + 2*lambda^T*cov(M_{v,x}^T,b_{v,x}) + V[b_{v,x}]
                % =E[M_{v,x}*lambda*lambda^T*M_{v,x}^T] - (E[M_{v,x}*lambda])^2 + 2*lambda^T*(E[M_{v,x}^T*b_{v,x}] - E[M_{v,x}^T]*E[b_{v,x}]) + V[b_{v,x}]
                V_Mvx_lambda_plus_bvx = E_Mvx_lambda_lambda_Mvx - (E_M_v_x*lambda_vec)^2 + 2*lambda_vec'*(E_Mvx_bvx - E_M_v_x'*E_b_Drx) + V_b_Drx;
                
                V_Phi(1) = lambda_tx^2*V_Mvx_lambda_plus_bvx;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                %dV_Phi(1) = ...
                
                % LCP variance matrix of V_Mvy_lambda_plus_bvy
                %fourth order expectation
                [E_ynyn,dE_ynyn_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_ynyx,dE_ynyx_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_ynyy,dE_ynyy_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                [E_yxyn,dE_yxyn_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_yxyx,dE_yxyx_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_yxyy,dE_yxyy_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                [E_yyyn,dE_yyyn_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,F,Fc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_yyyx,dE_yyyx_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,G,Gc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                [E_yyyy,dE_yyyy_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,H,Hc,Minv, Jg, dMdq, nq, Sigma_r, mu_r);
                
                % computing E_Mvy_lambda_lambda_Mvy
                E_Mvy_lambda_lambda_Mvy_quad = h^2*(E_ynyn*lambda_n^2 + E_ynyx*lambda_n*lambda_tx + E_ynyy*lambda_n*lambda_ty ...
                    + E_yxyn*lambda_tx*lambda_n + E_yxyx*lambda_tx^2 + E_yxyy*lambda_tx*lambda_ty ...
                    + E_yyyn*lambda_ty*lambda_n + E_yyyx*lambda_ty*lambda_tx + E_yyyy*lambda_ty^2);
                E_Mvy_lambda_lambda_Mvy_linear = 2*h*E_M_Dry_nr*lambda_n*gamma_left + 2*h*E_M_Dry_Drx*lambda_tx*gamma_left ...
                                          + 2*h*E_M_Dry_Dry*lambda_ty*gamma_left + gamma_left^2;
                E_Mvy_lambda_lambda_Mvy = E_Mvy_lambda_lambda_Mvy_quad + E_Mvy_lambda_lambda_Mvy_linear;
                
                % computing derivative of E_Mvy_lambda_lambda_Mvy
                dE_Mvy_lambda_lambda_Mvy_dh = 2*E_Mvy_lambda_lambda_Mvy_quad/h + (E_Mvy_lambda_lambda_Mvy_linear-gamma_left^2)/h;
                dE_Mvy_lambda_lambda_Mvy_dq0 = h^2/2*(dE_ynyn_dq*lambda_n^2 + dE_ynyx_dq*lambda_n*lambda_tx + dE_ynyy_dq*lambda_n*lambda_ty ...
                    + dE_yxyn_dq*lambda_tx*lambda_n + dE_yxyx_dq*lambda_tx^2 + dE_yxyy_dq*lambda_tx*lambda_ty ...
                    + dE_yyyn_dq*lambda_ty*lambda_n + dE_yyyx_dq*lambda_ty*lambda_tx + dE_yyyy_dq*lambda_ty^2) + h*dE_M_Dry_nr_dq*lambda_n*gamma_left ...
                    + h*dE_M_Dry_Drx_dq*lambda_tx*gamma_left + h*dE_M_Dry_Dry_dq*lambda_ty*gamma_left;
                dE_Mvy_lambda_lambda_Mvy_dv0 = zeros(nv,1);
                dE_Mvy_lambda_lambda_Mvy_dq1 = dE_Mvy_lambda_lambda_Mvy_dq0;
                dE_Mvy_lambda_lambda_Mvy_dv1 = zeros(nv,1);
                dE_Mvy_lambda_lambda_Mvy_du = zeros(nu,1);
                dE_Mvy_lambda_lambda_Mvy_dlambda_n = h^2*(2*E_ynyn*lambda_n + E_ynyx*lambda_tx + E_ynyy*lambda_ty + E_yxyn*lambda_tx + E_yyyn*lambda_ty) ...
                                                     + 2*h*E_M_Dry_nr*gamma_left;
                dE_Mvy_lambda_lambda_Mvy_dlambda_tx = h^2*(E_ynyx*lambda_n + E_yxyn*lambda_n + 2*E_yxyx*lambda_tx + E_yxyy*lambda_ty + E_yyyx*lambda_ty) ...
                                                      + 2*h*E_M_Dry_Drx*gamma_left;
                dE_Mvy_lambda_lambda_Mvy_dlambda_ty = h^2*(E_ynyy*lambda_n + E_yxyy*lambda_tx + E_yyyn*lambda_n + E_yyyx*lambda_tx + 2* E_yyyy*lambda_ty) ...
                                                      + 2*h*E_M_Dry_Dry*gamma_left;
                dE_Mvy_lambda_lambda_Mvy_dgamma = 2*h*E_M_Dry_nr*lambda_n + 2*h*E_M_Dry_Drx*lambda_tx ...
                                                    + 2*h*E_M_Dry_Dry*lambda_ty + 2*gamma_left;
                
                dE_Mvy_lambda_lambda_Mvy = [dE_Mvy_lambda_lambda_Mvy_dh;dE_Mvy_lambda_lambda_Mvy_dq0';dE_Mvy_lambda_lambda_Mvy_dv0;dE_Mvy_lambda_lambda_Mvy_dq1'; ...
                                            dE_Mvy_lambda_lambda_Mvy_dv1;dE_Mvy_lambda_lambda_Mvy_du;dE_Mvy_lambda_lambda_Mvy_dlambda_n;dE_Mvy_lambda_lambda_Mvy_dlambda_tx; ...
                                            dE_Mvy_lambda_lambda_Mvy_dgamma;zeros(3,1);dE_Mvy_lambda_lambda_Mvy_dgamma;0];
                 
                
                %Computing cov(M_{v,y}^T,b_{v,y}) = E[(M_{v,y}^T - E(M_{v,y}^T))(b_{v,y}-E(b_{v,y}))] = E[M_{v,y}^T*b_{v,y}] - E[M_{v,y}^T]*E[b_{v,y}];
                
                % E[M_{v,y}^T*b_{v,y}]
                [E_Mvy_bvy_Dry_nr_Dry,dE_Mvy_bvy_Dry_nr_Dry_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(H,Hc,F,Fc,H,Hc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                [E_Mvy_bvy_Dry_Drx_Dry,dE_Mvy_bvy_Dry_Drx_Dry_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(H,Hc,G,Gc,H,Hc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                [E_Mvy_bvy_Dry_Dry_Dry,dE_Mvy_bvy_Dry_Dry_Dry_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot] = expectation_third_order_multiply(H,Hc,H,Hc,H,Hc,Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r);
                E_Mvy_bvy = [E_Mvy_bvy_Dry_nr_Dry;E_Mvy_bvy_Dry_Drx_Dry;E_Mvy_bvy_Dry_Dry_Dry;E_b_Dry];
                
                %  V[M_{v,y}*lambda+b_{v,y}] = V[M_{v,y}*lambda] + cov(M_{v,y}*lambda,b_{v,y}) + cov(b_{v,y},M_{v,y}*lambda) + V[b_{v,y}]
                % = V[M_{v,y}*lambda] + 2*lambda^T*cov(M_{v,y}^T,b_{v,y}) + V[b_{v,y}]
                % =E[M_{v,y}*lambda*lambda^T*M_{v,y}^T] - (E[M_{v,y}*lambda])^2 + 2*lambda^T*cov(M_{v,y}^T,b_{v,y}) + V[b_{v,y}]
                % =E[M_{v,y}*lambda*lambda^T*M_{v,y}^T] - (E[M_{v,y}*lambda])^2 + 2*lambda^T*(E[M_{v,y}^T*b_{v,y}] - E[M_{v,y}^T]*E[b_{v,y}]) + V[b_{v,y}]
                V_Mvy_lambda_plus_bvy = E_Mvy_lambda_lambda_Mvy - (E_M_v_y*lambda_vec)^2 + 2*lambda_vec'*(E_Mvy_bvy - E_M_v_y'*E_b_Dry) + V_b_Dry;
                
                
                
                V_Phi(2) = lambda_ty^2*V_Mvy_lambda_plus_bvy;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                %dV_Phi(2) = ...
                
                function [E,dEdq,dEdqdot] = expectation_third_order_multiply(Ain, ain, Bin, bin, Cin, cin, Minv, Jg, J_blk, B, u, C, dMdq, dCdq, dCdqdot, nq, h, Sigma_r, mu_r)
                    %refactor inputs matrix
                    AAin = h*Minv*Jg'*Ain;
                    aain = h*Minv*Jg'*ain;
                    BBin = Jg'*Bin;
                    bbin = Jg'*bin;
                    CCin = J_blk'*Cin;
                    ccin = J_blk'*cin;
                    
                    % expectation
                    E = ((AAin*mu_r+aain)'*BBin+(BBin*mu_r+bbin)'*AAin)*Sigma_r*CCin'...
                        +(trace(AAin*Sigma_r*BBin')+(AAin*mu_r+aain)'*(BBin*mu_r+bbin))*(CCin*mu_r+ccin)';
                    
                    % derivative of expectation (d/dM * dM/dq)
                    % first term
                    term1 = ((AAin*mu_r+aain)'*BBin+(BBin*mu_r+bbin)'*AAin);
                    term2 = CCin';
                    dEdM = - Minv*(AAin*mu_r+aain)*term2'*Sigma_r*BBin'*Minv - Minv*(BBin*mu_r+bbin)*term2'*Sigma_r*AAin' - Minv*h*Jg'*Cin*Sigma_r*term1'*(B*u-C)'*Minv;
                    % second term
                    term3 = trace(AAin*Sigma_r*BBin')+(AAin*mu_r+aain)'*(BBin*mu_r+bbin);
                    term4 = (CCin*mu_r+ccin)';
                    dEdM = dEdM -Minv*term4*BBin*Sigma_r*AAin' -(AAin*mu_r+aain)*term4*(BBin*mu_r+bbin)'*Minv - Minv*(Jg'*Cin*mu_r+Jg'*cin)*term3*h*(B*u-C)'*Minv;
                    
                    dEdC = -(term1*Sigma_r*Cin'*Jg*h*Minv)';
                    dEdC = dEdC - (term3*(mu_r'*Cin'+cin')*Jg*Minv*h)';
                    
                    for i=1:nq
                        dEdq(i) = trace(dEdM'*dMdq(:,:,i)) + trace(dEdC'*dCdq(:,i));
                        dEdqdot(i) = trace(dEdC'*dCdqdot(:,i));
                    end
                    %[ToDo] take gradient w.r.t C and Jq
                end
                
                function [E,dEdq] = expectation_fourth_order_multiply(Ain, ain, Bin, bin, Cin, cin, Din, din, Minv, Jg, dMdq, nq, Sigma_r, mu_r)
                    %refactor inputs matrix
                    AAin = Minv*Jg'*Ain;
                    aain = Minv*Jg'*ain;
                    BBin = Jg'*Bin;
                    bbin = Jg'*bin;
                    CCin = Minv*Jg'*Cin;
                    ccin = Minv*Jg'*cin;
                    DDin = Jg'*Din;
                    ddin = Jg'*din;
                    
                    % expectation
                    E = trace(AAin*Sigma_r*(CCin'*DDin + DDin'*CCin)*Sigma_r*BBin') + ((AAin*mu_r+aain)'*BBin+(BBin*mu_r+bbin)'*AAin)*Sigma_r*(CCin'*(DDin*mu_r+ddin)+DDin'*(CCin*mu_r+ccin))...
                        +(trace(AAin*Sigma_r*BBin')+(AAin*mu_r+aain)'*(BBin*mu_r+bbin))*(trace(CCin*Sigma_r*DDin')+(CCin*mu_r+ccin)'*(DDin*mu_r+ddin));
                    
                    % derivative of expectation (d/dM * dM/dq)
                    % first term
                    dEdM = -Minv*BBin*Sigma_r*(CCin'*DDin+DDin'*CCin)*Sigma_r*AAin' - CCin*Sigma_r*AAin'*BBin*Sigma_r*DDin'*Minv - Minv*DDin*Sigma_r*AAin'*BBin*Sigma_r*CCin';
                    % second term
                    term1 = (CCin'*(DDin*mu_r+ddin)+DDin'*(CCin*mu_r+ccin));
                    term2 = ((AAin*mu_r+aain)'*BBin+(BBin*mu_r+bbin)'*AAin);
                    dEdM = dEdM - Minv*(AAin*mu_r+aain)*term1'*Sigma_r*BBin'*Minv - Minv*(BBin*mu_r+bbin)*term1'*Sigma_r*AAin' - CCin*Sigma_r*term2'*(DDin*mu_r+ddin)'*Minv ...
                         - Minv*DDin*Sigma_r*term2'*(mu_r'*CCin'+ccin');
                    % second term
                    term3 = trace(AAin*Sigma_r*BBin')+(AAin*mu_r+aain)'*(BBin*mu_r+bbin);
                    term4 = trace(CCin*Sigma_r*DDin')+(CCin*mu_r+ccin)'*(DDin*mu_r+ddin);
                    dEdM = dEdM -Minv*term4*BBin*Sigma_r*AAin' -(AAin*mu_r+aain)*term4*(BBin*mu_r+bbin)'*Minv - Minv*term3*DDin*Sigma_r*CCin'...
                         - (CCin*mu_r+ccin)*term3*(DDin*mu_r+ddin)'*Minv;
                    
                    %dEdC = 0;
                    
                    for i=1:nq
                        dEdq(i) = trace(dEdM'*dMdq(:,:,i));
                        dEdqot(i) = 0;
                    end
                    
                end
                
                
                
                % second foot, to be modified
                z_right_foot = [lambda(4:6);gamma(2)]';
                E_M_v_x = [h*E_M_Drx_nr, h*E_M_Drx_Drx, h*E_M_Drx_Dry, 1];
                E_M_v_y = [h*E_M_Dry_nr, h*E_M_Dry_Drx, h*E_M_Dry_Dry, 1];
                E_Phi(3) = lambda(2)*(E_M_v_x*z_right_foot + E_b_Drx);
                E_Phi(4) = lambda(3)*(E_M_v_y*z_right_foot + E_b_Dry);
                % to be modified
                
                g = obj.W*gamma + obj.M*lambda + obj.r;
                dg = [obj.M obj.W];
                
                % probabilistic version
                sigma = 0.5;
                lambda_parallel = [lambda(2);lambda(3);lambda(5);lambda(6)];
                
                %mu_cov = sigma^2*[lambda(1);lambda(4)];
                delta = 100;% coefficient
                f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);% - slack_var*ones(zdim,1);
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
                for j=1:obj.nD
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
            if nu>0
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