classdef RobustContactImplicitTrajectoryOptimization_Brick < DirectTrajectoryOptimization_Brick
    % phi, lambda
    properties
        nC
        nD % number of friction elements per contact
        nq
        nv
        nu
        nx
        nlambda
        lambda_inds
        cached_Px
        
        l_inds % orderered [lambda_N;lambda_f1;lambda_f2;...;gamma] for each contact sequentially
        lfi_inds % nD x nC indexes into lambda for each time step
        LCP_slack_inds % slack variable for LCP component <z,f(x,z) = LCP_slack > = 0
        lambda_mult
        ljl_inds  % joint limit forces
        jl_lb_ind  % joint indices where the lower bound is finite
        jl_ub_ind % joint indices where the lower bound is finite
        nJL % number of joint limits = length([jl_lb_ind;jl_ub_ind])
        F_ext_inds % index of external horizontal force
        nFext % number of external forces
        
        nonlincompl_constraints
        nonlincompl_slack_inds
        
        nonlincompl_constraints_purturb
        nonlincompl_slack_inds_purturb
        
        % complementarity matrix
        W
        r
        M
        C
        B
        u
        Jg
        Minv
        J_blk
        qdot_blk
        Sigma_r
        mu_r
        h
        
        dMdq
        dCdq
        dCdqdot
        dJgdq_i
        
        %debugging var
        verbose_print
    end
    
    properties (Constant)
        FORWARD_EULER = 1;
        BACKWARD_EULER = 2;
        MIDPOINT = 3;  % DEFAULT
        MIXED = 4;   % matched to TimeSteppingRigidBodyManipulator. Forward on qd, backward on q
    end
    
    methods
        function obj = RobustContactImplicitTrajectoryOptimization_Brick(plant,N,duration,options)
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
                options.integration_method = RobustContactImplicitTrajectoryOptimization_Brick.MIDPOINT;
            end
            
            obj = obj@DirectTrajectoryOptimization_Brick(plant,N,duration,options);
            
        end
        
        function obj = addDynamicConstraints(obj)
            nX = obj.plant.getNumStates();
            nU = obj.plant.getNumInputs();
            nq = obj.plant.getNumPositions();
            nL = obj.nC*(1+obj.nD);
            nFext = obj.nFext;
            N = obj.N;
            obj.nx = nX;
            obj.nu = nU;
            obj.nlambda = nL;
            
            constraints = cell(N-1,1);
            lincompl_constraints = cell(N-1,1);
            obj.nonlincompl_constraints = cell(N-1,1);
            obj.nonlincompl_slack_inds = cell(N-1,1);
            jlcompl_constraints = cell(N-1,1);
            dyn_inds = cell(N-1,1);
            
            n_vars = 2*nX + nU + 1 + obj.nC*(2+obj.nD) + obj.nJL + obj.nFext;
            cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.dynamics_constraint_fun);
            q0 = getZeroConfiguration(obj.plant);
            
            [~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(q0,false,obj.options.active_collision_options);
            
            external_force_max = 10;% random chosen value
            external_force_cnstr = BoundingBoxConstraint([-external_force_max],[external_force_max]);
            
            for i=1:obj.N-1,
                dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.l_inds(:,i);obj.ljl_inds(:,i);obj.F_ext_inds(:,i)};
                constraints{i} = cnstr.setName(sprintf('dynamics[%d]',i));
                
                obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                
                %constraint for bounded external force
                external_force_constraints{i} = external_force_cnstr.setName(sprintf('externalForce[%d]',i));
                obj = obj.addConstraint(external_force_constraints{i}, [obj.F_ext_inds(:,i)]);
                
                if obj.nC > 0
                    % indices for (i) gamma
                    gamma_inds = obj.l_inds(obj.nD+2:obj.nD+2:end,i);
                    % awkward way to pull out these indices, for (i) lambda_N and
                    % lambda_f
                    lambda_inds = obj.l_inds(repmat((1:1+obj.nD)',obj.nC,1) + kron((0:obj.nC-1)',(2+obj.nD)*ones(obj.nD+1,1)),i);
                    obj.lambda_inds(:,i) = lambda_inds;
                    
                    obj.options.nlcc_mode = 5;% robust mode
                    obj.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode);
                    obj.nonlincompl_slack_inds{i} = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraints{i}.n_slack;
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
                    
                    lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode);
                    obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds;obj.LCP_slack_inds(:,i)]);
                    
                    % add ERM cost for sliding velocity constraint uncertainty
                    %obj = obj.addCost(FunctionHandleObjective(2*nX+nU+6+2+1,@(h,x0,x1,u,lambda,gamma,verbose_print)ERMcost_slidingVelocity(obj,h,x0,x1,u,lambda,gamma),1), ...
                    %    {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);lambda_inds;gamma_inds});
                    
                    % add ERM cost for friction cone coefficient uncertainty
                    %obj = obj.addCost(FunctionHandleObjective(obj.nC*(2+obj.nD),@(lambda,gamma)ERMcost_friction(obj,lambda,gamma),1),{lambda_inds;gamma_inds});
                    
                    % (simple version) add ERM cost for normal distance uncertainty
                    %obj = obj.addCost(FunctionHandleObjective(nX+obj.nC*(1+obj.nD),@(x,lambda)ERMcost_normaldistance_Gaussian(obj,x,lambda),1), ...
                    %      {obj.x_inds(:,i+1);lambda_inds});
                end
                
                if obj.nJL > 0
                    % joint limit linear complementarity constraint
                    % lambda_jl /perp [q - lb_jl; -q + ub_jl]
                    W_jl = zeros(obj.nJL);
                    [r_jl,M_jl] = jointLimitConstraints(obj.plant,q0);
                    jlcompl_constraints{i} = LinearComplementarityConstraint_original(W_jl,r_jl,M_jl,obj.options.lincc_mode);
                    
                    obj = obj.addConstraint(jlcompl_constraints{i},[obj.x_inds(1:nq,i+1);obj.ljl_inds(:,i);obj.LCP_slack_inds(:,i)]);
                end
            end
            
            % robust variance cost with feedback control
            x_inds_stack = reshape(obj.x_inds,obj.N*nX,[]);
            u_inds_stack = reshape(obj.u_inds,obj.N*nU,[]);
            Fext_inds_stack = reshape(obj.F_ext_inds,obj.N*nFext,[]);
            lambda_inds_stack = reshape(obj.lambda_inds,(obj.N-1)*nL,[]);
            obj.cached_Px = zeros(obj.nx,obj.nx,obj.N);
            obj.cached_Px(:,:,1) = eye(obj.nx); %[ToDo: To be modified]
            
            obj = obj.addCost(FunctionHandleObjective(obj.N*(nX+nU),@(x_inds,Fext_inds)robustVariancecost(obj,x_inds,Fext_inds),1),{x_inds_stack;Fext_inds_stack});
            
            if (obj.nC > 0)
                obj = obj.addCost(FunctionHandleObjective(length(obj.LCP_slack_inds),@(slack)robustLCPcost(obj,slack),1),obj.LCP_slack_inds(:));
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
                
                f_numeric = f;
                df_numeric = df;
                %disp('check gradient')
                
                %                 if(f > 1e-3 || max(any(df > 1e-3)) == 1)
                %                     disp('come here')
                %                 end
                
                %                 [f_numeric,df_numeric] = geval(@(y) nonlincompl_fun_check(y),y,struct('grad_method','numerical'));
                %                 valuecheck(df,df_numeric,1e-3);
                %                 valuecheck(f,f_numeric,1e-3);
                
                function [f_num,df_num] = nonlincompl_fun_check(y)
                    nq = obj.plant.getNumPositions;
                    nv = obj.plant.getNumVelocities;
                    x = y(1:nq+nv+obj.nC);
                    z = y(nq+nv+obj.nC+1:end);
                    gamma = x(nq+nv+1:end);
                    
                    q = x(1:nq);
                    v = x(nq+1:nq+nv);
                    
                    [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    
                    f_num = zeros(obj.nC*(1+obj.nD),1);
                    df_num = zeros(obj.nC*(1+obj.nD),nq+nv+obj.nC*(2+obj.nD));
                    
                    f_num(1:1+obj.nD:end) = phi;
                    df_num(1:1+obj.nD:end,1:nq) = n;
                    for j=1:obj.nD,
                        f_num(1+j:1+obj.nD:end) = gamma+D{j}*v;
                        df_num(1+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                        df_num(1+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                        df_num(1+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
                    end
                end
            end
        end
        
        function [c,dc] = robustVariancecost(obj, x_full, Fext_full)
            
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nu = obj.nFext;%obj.plant.getNumInputs;
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px(:,:,1) = obj.cached_Px(:,:,1);
            
            % disturbance variance
            % currently only consider terrain height and friction coefficient
            Pw = diag([0.01, 0.04]); %[to be tuned]
            w_phi = normrnd(zeros(1,obj.N),sqrt(Pw(1,1)),1,obj.N);%height noise
            w_mu = normrnd(zeros(1,obj.N),sqrt(Pw(2,2)),1,obj.N);%friction coefficient noise
            w_noise = [w_phi;w_mu];
            
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            nw = size(Pw,1);
            K = [10*ones(nu,nq),ones(nu,nv)];
            
            %initialize c and dc
            kappa = 1;
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            c = kappa*trace(Px(:,:,1));
            dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
            
            % initialize gradient of Tr(V) w.r.t state vector x
            dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
            dTrVdu(:,:,1) = zeros(obj.N-1,nu);
            
            % time counter
            tStart = tic;
            
            for k = 1:obj.N-1%[Ye: double check the index]
                %Generate sigma points from Px(i+1)
                %[the sequential way to be modified]
                % currently, only use the initial variance matrix for the
                % propogation
                if k == 1
                    [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                    if d
                        diverge  = k;
                        return;
                    end
                    S = scale*S;
                    Sig(:,:,k) = [S -S];
                    for j = 1:(2*(obj.nx+nw))
                        Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                    end
                    w_averg = 1/(2*(obj.nx+nw));
                    x_mean(:,k) = zeros(obj.nx,1);
                    for j = 1:(2*(obj.nx+nw))
                        x_mean(:,k) = x_mean(:,k) + w_averg*Sig(1:obj.nx,j,k);
                    end
                    c = c + norm(x(:,k)-x_mean(:,k))^2;
                end
                
                %Propagate sigma points through nonlinear dynamics
                for j = 1:(2*(obj.nx+nw))
                    % a hacky way to implement the control input
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    Bmatrix(:,:,j,k) = [1;zeros(5,1)];%B;hand coding
                    
                    %                     [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j),false,obj.options.active_collision_options);
                    %                     J = zeros(nl,nq);
                    %                     J(1:1+obj.nD:end,:) = n;
                    %                     dJ = zeros(nl*nq,nq);
                    %                     dJ(1:1+obj.nD:end,:) = dn;%[double check how dn is factorized]
                    %
                    %                     for kk=1:length(D)
                    %                         J(1+kk:1+obj.nD:end,:) = D{kk};
                    %                         dJ(1+kk:1+obj.nD:end,:) = dD{kk};
                    %                     end
                    
                    % add feedback control
                    t = obj.plant.timestep*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                    [xdn,df] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                    
                    Sig(1:obj.nx/2,j,k+1) = xdn(1:obj.nx/2);
                    Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn(obj.nx/2+1:obj.nx);
                    dfdu(:,:,j,k+1) = [obj.plant.timestep^2*Hinv(:,:,j,k)*Bmatrix(:,:,j,k);obj.plant.timestep*Hinv(:,:,j,k)*Bmatrix(:,:,j,k)];
                    dfdSig(:,:,j,k+1) = df(:,2:obj.nx+1) - dfdu(:,:,j,k+1)*K;
                    dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                end
                
                %Calculate mean and variance w.r.t. [x_k] from sigma points
                w_averg = 1/(2*(obj.nx+nw));
                x_mean(:,k+1) = zeros(obj.nx,1);
                for j = 1:(2*(obj.nx+nw))
                    x_mean(:,k+1) = x_mean(:,k+1) + w_averg*Sig(1:obj.nx,j,k+1);
                end
                Px(:,:,k+1) = zeros(obj.nx);
                %alpha = 1e-3;
                %w_coeff = (1/(2*alpha^2*(obj.nx+nw)));
                
                for j = 1:(2*(obj.nx+nw))
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
                end
                
                % accumulate returned cost
                c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;
                for j = 1:(2*(obj.nx+nw))
                    V_comp = (Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
                    c = c + kappa*trace(w*V_comp);
                end
                
                % derivative of variance matrix
                % gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,k+1) = zeros(obj.N-1,obj.nx);
                dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                for j=k:-1:1
                    dTrVdx(j,:,k+1) = zeros(1,obj.nx);
                    dTrVu(j,:,k+1) = zeros(1,nu);
                    dmeanRdx(j,:,k+1) = zeros(1,obj.nx);
                    dmeanRdu(j,:,k+1) = zeros(1,nu);
                    
                    % gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(obj.nx);
                    % gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(obj.nx,1);
                    
                    for i=1:2*(obj.nx+nw)
                        if i == 1
                            for m = 1:(2*(obj.nx+nw))% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                % gradient of Tr(V_{k+1}) w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(obj.nx);
                                dSig_m_kplus1_du = zeros(obj.nx,1);
                                
                                chain_rule_indx = k-j;
                                if j ~= 1
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2)*eye(obj.nx);
                                end
                                dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                
                                while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                    chain_rule_indx = chain_rule_indx - 1;
                                end
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run 2*(obj.nx+nw) times in total to obtain
                        % gradient w.r.t sigma points
                        dSig_i_kplus1_dx = zeros(obj.nx);
                        dSig_i_kplus1_du = zeros(obj.nx,1);
                        chain_rule_indx = k-j;
                        if j ~= 1
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                        else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2)*eye(obj.nx);
                        end
                        dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                            dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                        dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_averg*dSig_m_kplus1_dx_sum);
                        dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_averg*dSig_m_kplus1_du_sum);
                    end
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_dx_sum);
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_du_sum);
                end
            end
            tElapsed = toc(tStart);
            
            dc = [];
            % cost gradient w.r.t x at first time step is zero
            for k=1:obj.N % index for x_k
                dTrV_sum_dx_k = zeros(1, obj.nx);
                if (k == 1)
                    dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                else
                    dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';
                end
                for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(k,:,kk);
                    dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(k,:,kk);
                end
                dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            % c(k=1) = Px(:,:,1)
            for k=1:obj.N % index for u_k
                dTrV_sum_du_k = zeros(1, nu);
                dmeanR_sum_du_k = zeros(1, nu);
                for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(k,:,kk);
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(k,:,kk);
                end
                dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
            end
            
            obj.cached_Px = Px;
            
            % check gradient
            disp('check gradient')
            c_numeric = c;
            dc_numeric = dc;
            
            X0 = [x_full; Fext_full];
            %X0 = X0 + randn(size(X0))*0.1;
            
            fun = @(X0) robustVariancecost_check(obj, X0)
            DerivCheck(fun, X0)
            
            u0 = Fext_full;
            fun = @(u) robustVariancecost_check_uonly(obj, u, X0)
            DerivCheck(fun, u0)
            
            x0 = x_full;
            fun = @(x) robustVariancecost_check_xonly(obj, x, X0)
            DerivCheck(fun, x0)
            
            [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_check(obj,X0),X0,struct('grad_method','numerical'));
            valuecheck(dc,dc_numeric,1e-5);
            valuecheck(c,c_numeric,1e-5);
            
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
                
                % Pick a random small vector in parameter space
                tol = 1e-6;  % Size of numerical step to take
                rr = sqrt(eps(X0));%randn(length(X0),1)*tol;  % Generate small random-direction vector
                
                % Evaluate at symmetric points around X0
                f1 = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0
                f2 = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0
                
                % Print results
                fprintf('Derivs: Analytic vs. Finite Diff = [%.12e, %.12e]\n', dot(rr, JJ), f2-f1);
                dd =  dot(rr, JJ)-f2+f1
            end
            
            function [c,dc] = robustVariancecost_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nu = obj.nFext;%obj.plant.getNumInputs;
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px(:,:,1) = obj.cached_Px(:,:,1);
                
                % disturbance variance
                % currently only consider terrain height and friction coefficient
                Pw = diag([0.01, 0.04]); %[to be tuned]
                w_phi = normrnd(zeros(1,obj.N),sqrt(Pw(1,1)),1,obj.N);%height noise
                w_mu = normrnd(zeros(1,obj.N),sqrt(Pw(2,2)),1,obj.N);%friction coefficient noise
                w_noise = [w_phi;w_mu];
                
                scale = .1;% [to be tuned]
                nw = size(Pw,1);
                K = [10*ones(nu,nq),ones(nu,nv)];
                
                %initialize c and dc
                kappa = 1;
                x_mean = zeros(obj.nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
                %c = kappa*trace(Px(:,:,1));
                c = 0;
                dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
                
                % initialize gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
                dTrVdu(:,:,1) = zeros(obj.N-1,nu);
                
                tStart = tic;
                
                for k = 1:obj.N-1%[Ye: double check the index]
                    %Generate sigma points from Px(i+1)
                    %[the sequential way to be modified]
                    % currently, only use the initial variance matrix for the
                    % propogation
                    if k == 1
                        [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                        if d
                            diverge  = k;
                            return;
                        end
                        S = scale*S;
                        Sig(:,:,k) = [S -S];
                        for j = 1:(2*(obj.nx+nw))
                            Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                        end
                        w_averg = 1/(2*(obj.nx+nw));
                        x_mean(:,k) = zeros(obj.nx,1);
                        for j = 1:(2*(obj.nx+nw))
                            x_mean(:,k) = x_mean(:,k) + w_averg*Sig(1:obj.nx,j,k);
                        end
                        c = c + norm(x(:,k)-x_mean(:,k))^2;
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:(2*(obj.nx+nw))
                        % a hacky way to implement the control input (only apply to fallingbrick example)
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        Bmatrix(:,:,j,k) = [1;zeros(5,1)];%B;hand coding
                        
                        % add feedback control
                        t = obj.plant.timestep*(k-1);
                        u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                        [xdn,df] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                        
                        Sig(1:obj.nx/2,j,k+1) = xdn(1:obj.nx/2);
                        Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn(obj.nx/2+1:obj.nx);
                        dfdu(:,:,j,k+1) = [obj.plant.timestep^2*Hinv(:,:,j,k)*Bmatrix(:,:,j,k);obj.plant.timestep*Hinv(:,:,j,k)*Bmatrix(:,:,j,k)];
                        dfdSig(:,:,j,k+1) = df(:,2:obj.nx+1) - dfdu(:,:,j,k+1)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                    end
                    
                    %Calculate mean and variance w.r.t. [x_k] from sigma points
                    w_averg = 1/(2*(obj.nx+nw));
                    x_mean(:,k+1) = zeros(obj.nx,1);
                    for j = 1:(2*(obj.nx+nw))
                        x_mean(:,k+1) = x_mean(:,k+1) + w_averg*Sig(1:obj.nx,j,k+1);
                    end
                    Px(:,:,k+1) = zeros(obj.nx);
                    %alpha = 1e-3;
                    %w_coeff = (1/(2*alpha^2*(obj.nx+nw)));
                    w = 0.5/scale^2;
                    for j = 1:(2*(obj.nx+nw))
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
                    end
                    
                    % accumulate returned cost
                    c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;
%                     for j = 1:(2*(obj.nx+nw))
%                         V_comp = (Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
%                         c = c + kappa*trace(w*V_comp);
%                     end
                    
                    % derivative of variance matrix
                    % gradient of Tr(V) w.r.t state vector x
                    dTrVdx(:,:,k+1) = zeros(obj.N-1,obj.nx);
                    dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                    dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                    dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                    
                    for j=k:-1:1
                        dTrVdx(j,:,k+1) = zeros(1,obj.nx);
                        dTrVu(j,:,k+1) = zeros(1,nu);
                        dmeanRdx(j,:,k+1) = zeros(1,obj.nx);
                        dmeanRdu(j,:,k+1) = zeros(1,nu);
                        
                        % gradient w.r.t state x
                        dSig_m_kplus1_dx_sum = zeros(obj.nx);
                        % gradient w.r.t control u
                        dSig_m_kplus1_du_sum = zeros(obj.nx,1);
                        
                        for i=1:2*(obj.nx+nw)
                            if i == 1
                                for m = 1:(2*(obj.nx+nw))% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                    % gradient of Tr(V_{k+1}) w.r.t control x and u
                                    dSig_m_kplus1_dx = zeros(obj.nx);
                                    dSig_m_kplus1_du = zeros(obj.nx,1);
                                    
                                    chain_rule_indx = k-j;
                                    if j ~= 1
                                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                    else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2)*eye(obj.nx);
                                    end
                                    dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                    
                                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                        dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                        dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                        chain_rule_indx = chain_rule_indx - 1;
                                    end
                                    dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                                end
                            end
                            
                            % run 2*(obj.nx+nw) times in total to obtain
                            % gradient w.r.t sigma points
                            dSig_i_kplus1_dx = zeros(obj.nx);
                            dSig_i_kplus1_du = zeros(obj.nx,1);
                            chain_rule_indx = k-j;
                            %if j ~= 1
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                            %else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            %    dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2)*eye(obj.nx);
                            %end
                            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            
                            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_averg*dSig_m_kplus1_dx_sum);
                            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_averg*dSig_m_kplus1_du_sum);
                        end
                        
                        % gradient of mean residual w.r.t state x and control u, assume norm 2
                        dmeanRdx(j,:,k+1) = 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_dx_sum);
                        dmeanRdu(j,:,k+1) = 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_du_sum);
                    end
                end
                tElapsed = toc(tStart);
                
                dc = [];
                % cost gradient w.r.t x at first time step is zero
                for k=1:obj.N % index for x_k
                    dTrV_sum_dx_k = zeros(1, obj.nx);
                    if (k == 1)
                        dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                    else
                        dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';
                    end
                    for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(k,:,kk);
                        dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(k,:,kk);
                    end
                    dc = [dc, dmeanR_sum_dx_k];%+kappa*dTrV_sum_dx_k
                end
                
                % cost gradient w.r.t u at first time step is zero, since
                % c(k=1) = Px(:,:,1)
                for k=1:obj.N % index for u_k
                    dTrV_sum_du_k = zeros(1, nu);
                    dmeanR_sum_du_k = zeros(1, nu);
                    for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(k,:,kk);
                        dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(k,:,kk);
                    end
                    dc = [dc, dmeanR_sum_du_k];%+kappa*dTrV_sum_du_k
                end
                
                obj.cached_Px = Px;
            end
        end
        
        function [c,dc] = robustLCPcost(obj, slack_var)
            c = 1000*sum(slack_var);
            dc = 1000*ones(1,length(slack_var));
        end
        
        function [f,df] = ERMcost_normaldistance_Gaussian(obj, x, lambda)
            
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            z = lambda;
            q = x(1:nq);
            v = x(nq+1:nq+nv);
            
            %                 % von Mises-Fisher distribution for quaternion rotation vector
            %                 mu_dirc = [1,0,0,0]'; % this is for flat terrain, no rotation. can be terrain dependent.
            %
            %                 kappa = 10;
            %                 I_kappa_plus = exp(kappa) + exp(-kappa);
            %                 I_kappa_minus = exp(kappa) - exp(-kappa);
            %
            %                 h_kappa = (kappa*I_kappa_plus - I_kappa_minus)/(kappa*I_kappa_minus);
            %                 mu_r = mu_dirc*h_kappa;
            %                 Sigma_r = h_kappa/kappa * eye(4) + (1 - 3*h_kappa/kappa - h_kappa^2)*mu_r*mu_r';
            %
            %                 %remove distribution and make it deterministic
            %                 mu_r=[1;0;0;0];
            %                 Sigma_r = zeros(4);
            %
            %                 obj.mu_r = mu_r;
            %                 obj.Sigma_r = Sigma_r;
            %
            %                 mu_w = mu_r(1);mu_x = mu_r(2);mu_y = mu_r(3);mu_z = mu_r(4);
            %
            %                 Rbar = [1-2*mu_y^2-2*mu_z^2, 2*mu_x*mu_y-2*mu_z*mu_w, 2*mu_x*mu_z+2*mu_y*mu_w;
            %                         2*mu_x*mu_y+2*mu_z*mu_w, 1-2*mu_x^2-2*mu_z^2, 2*mu_y*mu_z-2*mu_x*mu_w;
            %                         2*mu_x*mu_z-2*mu_y*mu_w, 2*mu_y*mu_z+2*mu_x*mu_w, 1-2*mu_x^2-2*mu_y^2];
            %
            %                 % normal direction n
            %                 F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
            %                     -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
            %                      2*mu_w, -2*mu_x, -2*mu_y, 2*mu_z];
            %                 Fc = Rbar(:,3) - F*mu_r;
            
            [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
            
            %                 phi_offset = n*q - phi;
            %                 f_vec = zeros(obj.nC,1);
            %                 df_vec = zeros(obj.nC,nq+nv+obj.nC*(2+obj.nD));
            
            %                 for contact_index = 1:obj.nC
            %                     Jg(1:obj.nD,:) = [D{1}(contact_index,:); D{2}(contact_index,:); D{3}(contact_index,:); D{4}(contact_index,:)];
            %                     Jg(1+obj.nD,:) = n(contact_index,:);
            %                 end
            
            %                 sigma_height = 0.05;
            %                 kappa = 1;% covariance coefficient
            %                 E_phi = (mu_r'*F' + Fc')*n*q + phi_offset;
            %                 V_phi = F*Sigma_r*F'*n*q;
            %                 dE_phi(:,1:nq) = (mu_r'*F' + Fc')*n;
            %                 dV_phi(:,1:nq) = F*Sigma_r*F'*n;
            
            
            lambda_n = z(1:1+obj.nD:end);
            
            % normal distance uncertainty
            mu_phi = zeros(obj.nC,1);
            Sigma_phi = 0.1*ones(obj.nC,1);
            kappa = 1;% covariance coefficient
            E_Phi = phi + mu_phi;
            V_Phi = Sigma_phi;
            E_Phi = E_Phi.*lambda_n;
            V_Phi = kappa*V_Phi.*lambda_n;
            
            dE_Phi = zeros(obj.nC,nq+nv+obj.nC*(1+obj.nD));
            dV_Phi = zeros(obj.nC,nq+nv+obj.nC*(1+obj.nD));
            
            for contact_index = 1:obj.nC
                dE_Phi(contact_index, 1:nq) = n(contact_index,:);
                dE_Phi(contact_index, nq+nv+1+(1+obj.nD)*(contact_index-1)) = phi(contact_index) + mu_phi(contact_index);
                dV_Phi(contact_index, nq+nv+1+(1+obj.nD)*(contact_index-1)) = kappa*Sigma_phi(contact_index);
            end
            
            delta = 1000;% coefficient
            f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);
            df = delta*(E_Phi'*dE_Phi + V_Phi'*dV_Phi);
            
            %                 persistent normal_LCP_robust_NCP
            %                 normal_LCP_robust_NCP = [normal_LCP_robust_NCP, f_vec.*[z(2);z(3);z(5);z(6)]];
            %                 if length(normal_LCP_robust_NCP) == obj.N-1
            %                     normal_LCP_robust_NCP
            %                     normal_LCP_robust_NCP = [];
            %                 end
            
            %disp('check gradient')
            f_numeric = f;
            df_numeric = df;
            
            %             if (f > 1e-3 || max(any(df > 1e-3)) == 1)
            %                 disp('come here')
            %             end
            
            %             [f_numeric,df_numeric] = geval(@(x,lambda) ERMcost_normaldistance_Gaussian_check(obj,x,lambda),x,lambda,struct('grad_method','numerical'));
            %             valuecheck(df,df_numeric,1e-3);
            %             valuecheck(f,f_numeric,1e-3);
            
            function [f,df] = ERMcost_normaldistance_Gaussian_check(obj, x, lambda)
                
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                z = lambda;
                q = x(1:nq);
                v = x(nq+1:nq+nv);
                
                [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                
                lambda_n = z(1:1+obj.nD:end);
                
                % normal distance uncertainty
                mu_phi = zeros(obj.nC,1);
                Sigma_phi = 0.1*ones(obj.nC,1);
                kappa = 1;% covariance coefficient
                E_Phi = phi + mu_phi;
                V_Phi = Sigma_phi;
                E_Phi = E_Phi.*lambda_n;
                V_Phi = kappa*V_Phi.*lambda_n;
                
                dE_Phi = zeros(obj.nC,nq+nv+obj.nC*(1+obj.nD));
                dV_Phi = zeros(obj.nC,nq+nv+obj.nC*(1+obj.nD));
                
                for contact_index = 1:obj.nC
                    dE_Phi(contact_index, nq+nv+1+(1+obj.nD)*(contact_index-1)) = phi(contact_index) + mu_phi(contact_index);
                    dV_Phi(contact_index, nq+nv+1+(1+obj.nD)*(contact_index-1)) = kappa*Sigma_phi(contact_index);
                end
                
                delta = 1000;% coefficient
                f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);
                df = delta*(E_Phi'*dE_Phi + V_Phi'*dV_Phi);
                
            end
        end
        
        function [f,df] = ERMcost_friction(obj, lambda, gamma)
            
            zdim = size(obj.W,2);
            xdim = size(obj.M,2);
            
            %  Wz+Mx+r = mu*lambda_N - sum(lambda_fi)
            g = obj.W*gamma + obj.M*lambda + obj.r;
            dg = [obj.M obj.W];
            
            %                 % distribution-free version
            %                 delta = 1000;% coefficient
            %                 f = delta/2 * norm(gamma.*g)^2;% - slack_var*ones(zdim,1);
            %                 df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
            
            % probabilistic version
            sigma = 0.5;
            lambda_n = lambda(1:1+obj.nD:end);
            mu_cov = sigma^2*lambda_n.^2;% make sure it is squared
            mu_mu = 1;
            Sigma_mu = 0.25;
            E_Phi = g.*gamma;
            V_Phi = Sigma_mu*lambda_n.*gamma;
            
            dE_Phi = zeros(obj.nC,obj.nC*(2+obj.nD));
            dV_Phi = zeros(obj.nC,obj.nC*(2+obj.nD));
            
            dE_Phi = [diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
            
            for contact_index = 1:obj.nC
                dV_Phi(contact_index, 1+(1+obj.nD)*(contact_index-1)) = Sigma_mu*gamma(contact_index);
                dV_Phi(contact_index, (1+obj.nD)*obj.nC+contact_index) = Sigma_mu*lambda_n(contact_index);
            end
            
            delta = 10^7;% penalty coefficient
            f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);
            df = delta*(E_Phi'*dE_Phi + V_Phi'*dV_Phi);
            
            
            %disp('check gradient')
            f_numeric = f;
            df_numeric = df;
            
            %             if (f > 1e-3 || max(any(df > 1e-3)) == 1)
            %                 disp('come here')
            %             end
            
            %             [f_numeric,df_numeric] = geval(@(lambda,gamma) ERMcost_friction_check(obj,lambda,gamma),lambda,gamma,struct('grad_method','numerical'));
            %             valuecheck(df,df_numeric,1e-15);
            %             valuecheck(f,f_numeric,1e-15);
            
            function [f,df] = ERMcost_friction_check(obj, lambda, gamma)
                
                zdim = size(obj.W,2);
                xdim = size(obj.M,2);
                
                %  Wz+Mx+r = mu*lambda_N - sum(lambda_fi)
                g = obj.W*gamma + obj.M*lambda + obj.r;
                dg = [obj.M obj.W];
                
                %                 % distribution-free version
                %                 delta = 1000;% coefficient
                %                 f = delta/2 * norm(gamma.*g)^2;% - slack_var*ones(zdim,1);
                %                 df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
                % probabilistic version
                sigma = 0.5;
                lambda_n = lambda(1:1+obj.nD:end);
                mu_cov = sigma^2*lambda_n.^2;% make sure it is squared
                mu_mu = 1;
                Sigma_mu = 0.25;
                E_Phi = g.*gamma;
                V_Phi = Sigma_mu*lambda_n.*gamma;
                
                dE_Phi = zeros(obj.nC,obj.nC*(2+obj.nD));
                dV_Phi = zeros(obj.nC,obj.nC*(2+obj.nD));
                
                dE_Phi = [diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
                
                for contact_index = 1:obj.nC
                    dV_Phi(contact_index, 1+(1+obj.nD)*(contact_index-1)) = Sigma_mu*gamma(contact_index);
                    dV_Phi(contact_index, (1+obj.nD)*obj.nC+contact_index) = Sigma_mu*lambda_n(contact_index);
                end
                
                delta = 10^7;% penalty coefficient
                f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);
                df = delta*(E_Phi'*dE_Phi + V_Phi'*dV_Phi);
            end
        end
        
        function [f,df] = ERMcost_slidingVelocity(obj, h, x0, x1, u, lambda, gamma)
            
            zdim = size(obj.W,2);
            xdim = size(obj.M,2);
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nu = obj.plant.getNumInputs;
            nl = length(lambda);
            
            obj.nq = nq;
            obj.nv = nv;
            obj.nu = nu;
            
            obj.h = h;
            
            % von Mises-Fisher distribution for quaternion rotation vector
            mu_dirc = [1,0,0,0]'; % this is for flat terrain, no rotation. can be terrain dependent.
            %mu_dirc = [0.2,0.15,0.3,0.9206]'; % for a tilt terrain
            
            kappa = 10;
            I_kappa_plus = exp(kappa) + exp(-kappa);
            I_kappa_minus = exp(kappa) - exp(-kappa);
            
            h_kappa = (kappa*I_kappa_plus - I_kappa_minus)/(kappa*I_kappa_minus);
            mu_r = mu_dirc*h_kappa;
            Sigma_r = h_kappa/kappa * eye(4) + (1 - 3*h_kappa/kappa - h_kappa^2)*mu_r*mu_r';
            
            %remove distribution and make it deterministic
            mu_r=[1;0;0;0];
            Sigma_r = zeros(4);
            
            obj.mu_r = mu_r;
            obj.Sigma_r = Sigma_r;
            
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
            obj.C = C;
            obj.B = B;
            obj.u = u;
            obj.Minv = Minv;
            
            for i=1:nq
                dMdq(:,:,i) = reshape(dM(:,i),[nq, nq]);
                dCdq(:,i) = reshape(dC(:,i),[nq, 1]);
                dCdqdot(:,i) = reshape(dC(:,i+nq),[nq, 1]);
            end
            
            obj.dMdq = dMdq;
            obj.dCdq = dCdq;
            obj.dCdqdot = dCdqdot;
            
            qdot_prev = v0;
            u_prev = u;
            
            if nl>0
                [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                % construct J and dJ from n,D,dn, and dD so they relate to the
                % lambda vector
                J = zeros(nl,nq);
                J(1:1+obj.nD:end,:) = n;
                dJ = zeros(nl*nq,nq);
                dJ(1:1+obj.nD:end,:) = dn;%[double check how dn is factorized]
                
                for j=1:length(D),
                    J(1+j:1+obj.nD:end,:) = D{j};
                    dJ(1+j:1+obj.nD:end,:) = dD{j};
                end
            end
            
            E_Phi = zeros(4,1);
            V_Phi = zeros(4,1);
            
            for foot_indx = 1:2
                Jg(1:2,:) = J(3*foot_indx-1:3*foot_indx,:);
                Jg(3,:) = J(3*foot_indx-2,:);
                
                obj.Jg = Jg;
                
                for i=1:nq
                    for j =1:((1+obj.nD)*2)
                        if foot_indx == 1
                            dJgdq(:,j,i) = dJ(2*(1+obj.nD)*(j-1)+1:2*(1+obj.nD)*(j-1)+3,i);
                        elseif foot_indx == 2
                            dJgdq(:,j,i) = dJ(2*(1+obj.nD)*(j-1)+4:2*(1+obj.nD)*(j-1)+6,i);
                        end
                    end
                end
                
                % Composed matrices
                Fy = H;
                Fyc = Hc;
                U = Minv*Jg'*F;
                V = U'*Jg'*G;
                Vy = U'*Jg'*Fy;
                Z = Minv*Jg'*Fc;
                Zy = Minv*Jg'*Fyc;
                X = G'*Jg*Z;
                Xy = Fy'*Jg*Z;
                Xxy = G'*Jg*Zy;
                K = Minv*Jg'*G;
                Ky = Minv*Jg'*Fy;
                L = Minv*Jg'*Gc;
                Ly = Minv*Jg'*Fyc;
                Wx = K'*Jg'*G;
                Wy = Ky'*Jg'*Fy;
                Wxy = K'*Jg'*Fy;
                Y = U'*Jg'*Fc;
                Yx = K'*Jg'*Gc;
                Yy = Ky'*Jg'*Fyc;%[double check]
                Q = U'*Jg'*Gc;
                Qy = U'*Jg'*Fyc;
                Qxy = Ky'*Jg'*Gc;%[double check]
                O = (V+V')*mu_r+X+Q;
                Oy = (Vy+Vy')*mu_r+Xy+Qy;
                Oxy = (Wxy+Wxy')*mu_r+Xxy+Qxy;
                
                qdot_blk = (qdot_prev + Minv*(B*u_prev - C)*h);
                J_blk = Jg*qdot_blk;
                T = G'*J_blk;
                Tc = Gc'*J_blk;
                Ty = Fy'*J_blk;
                Tyc = Fyc'*J_blk;
                obj.qdot_blk = qdot_blk;
                obj.J_blk = J_blk;
                
                %--------------- second LCP condition ---------------%
                % expectation and covariance of M_v_x
                E_M_Drx_nr = trace(V*Sigma_r) + mu_r'*V*mu_r + Z'*Jg'*(G*mu_r + Gc) + mu_r'*U'*Jg'*Gc;
                
                V_M_Drx_nr = trace(U*Sigma_r*(V+V')*Sigma_r*G'*Jg) + O'*Sigma_r*O ...
                    +(trace(U*Sigma_r*G'*Jg)+mu_r'*V*mu_r+mu_r'*Q+X'*mu_r+Z'*Jg'*Gc)^2 - E_M_Drx_nr^2;
                
                E_M_Drx_Drx = trace(Wx*Sigma_r) + mu_r'*Wx*mu_r + L'*Jg'*(2*G*mu_r + Gc);
                V_M_Drx_Drx = 2*trace(K*Sigma_r*Wx*Sigma_r*G'*Jg) + 4*(mu_r'*Wx + Yx')*Sigma_r*(Wx*mu_r + Yx) ...
                    +(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc)^2 - E_M_Drx_Drx^2; %[not used]
                
                E_M_Drx_Dry = trace(Wxy*Sigma_r) + mu_r'*Wxy*mu_r + L'*Jg'*(Fy*mu_r + Fyc) + mu_r'*K'*Jg'*Fyc;
                E_M_Dry_Drx = E_M_Drx_Dry;
                % no V_M_Drx_Dry defined, not used
                
                for i=1:nq
                    % expectation derivative w.r.t q and qdot
                    obj.dJgdq_i = dJgdq(:,:,i);
                    dE_M_Drx_nr_dq(i,:) = trace( (-Minv*Jg'*(F*(Sigma_r + mu_r*mu_r')*G' + Fc*(G*mu_r + Gc)' + F*mu_r*Gc')*Jg*Minv)'*dMdq(:,:,i) ) ...
                        + jacobian_gradient(F',Minv,G*Sigma_r) + jacobian_gradient(mu_r'*F',Minv,G*mu_r) ...
                        + jacobian_gradient(Fc',Minv,G*mu_r + Gc) + jacobian_gradient(mu_r'*F',Minv,Gc);
                    dE_M_Drx_nr_dqdot(i,:) = 0;
                    
                    dE_M_Drx_Drx_dq(i,:) = trace( (-Minv*Jg'*(G*(Sigma_r + mu_r*mu_r')*G' + Gc*(2*G*mu_r + Gc)')*Jg*Minv)'*dMdq(:,:,i) ) ...
                        + jacobian_gradient(G',Minv,G*Sigma_r) + jacobian_gradient(mu_r'*G',Minv,G*mu_r) ...
                        + jacobian_gradient(Gc',Minv,2*G*mu_r+Gc);
                    dE_M_Drx_Drx_dqdot(i,:) = 0;
                    
                    dE_M_Drx_Dry_dq(i,:) = trace( (-Minv*Jg'*(Fy*(Sigma_r + mu_r*mu_r')*G' + Fyc*(G*mu_r + Gc)' + Fy*mu_r*Gc')*Jg*Minv)'*dMdq(:,:,i) ) ...
                        + jacobian_gradient(Fy',Minv,G*Sigma_r) + jacobian_gradient(mu_r'*Fy',Minv,G*mu_r) ...
                        + jacobian_gradient(Fyc',Minv,G*mu_r + Gc) + jacobian_gradient(mu_r'*Fy',Minv,Gc);
                    dE_M_Drx_Dry_dqdot(i,:) = 0;
                    
                    dE_M_Dry_Drx_dq(i,:) = dE_M_Drx_Dry_dq(i,:);
                    dE_M_Dry_Drx_dqdot(i,:) = 0;
                    
                    % covariance derivative w.r.t q and qdot
                    dV_M_Drx_nr_dq_first_chain(:,:,i) = -K*Sigma_r*(V+V')*Sigma_r*U' - U*Sigma_r*V*Sigma_r*K' - K*Sigma_r*V*Sigma_r*U' ...
                        -U*(mu_r*O'+O*mu_r')*K'-K*(mu_r*O'+O*mu_r')*U'-Z*O'*K'-L*O'*U'-K*O*Z'-U*O*L' ...
                        +2*(trace(U*Sigma_r*G'*Jg)+mu_r'*V*mu_r+mu_r'*Q+X'*mu_r+Z'*Jg'*Gc) ...
                        *(-K*Sigma_r*U' - U*mu_r*mu_r'*K' - U*mu_r*L' - Z*mu_r'*K' - Z*L');
                    dV_M_Drx_nr_dq_first_chain(:,:,i) = dV_M_Drx_nr_dq_first_chain(:,:,i)';% transponse due to that the equation above is derived based on dV_M_nr_Drx
                    dV_M_Drx_nr_dq(i,:) = trace(dV_M_Drx_nr_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_nr_dq(i,:) = dV_M_Drx_nr_dq(i,:) - 2*E_M_Drx_nr*dE_M_Drx_nr_dq(i,:);%[double check this part]
                    dV_M_Drx_nr_dq(i,:) = dV_M_Drx_nr_dq(i,:) + jacobian_gradient3(Minv,F*Sigma_r*F',Minv,G*Sigma_r*G',eye(nq)) ...
                        + jacobian_gradient2(mu_r'*F'+Fc',Minv,G*Sigma_r*F',Minv,G*mu_r+Gc) ...
                        + jacobian_gradient2(mu_r'*F'+Fc',Minv,G*Sigma_r*G',Minv,F*mu_r+Fc) ...
                        + jacobian_gradient2(mu_r'*G'+Gc',Minv,F*Sigma_r*F',Minv,G*mu_r+Gc) ...
                        + jacobian_gradient2(mu_r'*G'+Gc',Minv,F*Sigma_r*G',Minv,F*mu_r+Fc) ...
                        +2*(trace(U*Sigma_r*G'*Jg)+mu_r'*V*mu_r+mu_r'*Q+X'*mu_r+Z'*Jg'*Gc) ...
                        *(jacobian_gradient1(Minv,F*Sigma_r*G') + jacobian_gradient(mu_r'*F',Minv,G*mu_r) ...
                        + jacobian_gradient(mu_r'*F',Minv,Gc) + jacobian_gradient(Fc',Minv,G*mu_r+Gc));
                    
                    dV_M_Drx_nr_dqdot(i,:) = 0;
                    
                    dV_M_Drx_Drx_dq_first_chain(:,:,i) = -4*K*Sigma_r*Wx*Sigma_r*K' + 4*(-K*mu_r*mu_r'*Wx*Sigma_r*K' - K*Sigma_r*Wx*mu_r*mu_r'*K' ...
                        -K*mu_r*Yx'*Sigma_r*K'-K*Sigma_r*Wx*mu_r*L'-L*mu_r'*Wx*Sigma_r*K'-K*Sigma_r*Yx*mu_r'*K'-L*Yx'*Sigma_r*K'-K*Sigma_r*Yx*L') ...
                        +2*(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc) ...
                        *(-K*(Sigma_r + mu_r*mu_r')*K' - 2*K*mu_r*L' - L*L');
                    dV_M_Drx_Drx_dq(i,:) = trace(dV_M_Drx_Drx_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_Drx_dq(i,:) = dV_M_Drx_Drx_dq(i,:) - 2*E_M_Drx_Drx*dE_M_Drx_Drx_dq(i,:);%[double check this part]
                    dV_M_Drx_Drx_dq(i,:) = dV_M_Drx_Drx_dq(i,:) + 2*jacobian_gradient3(Minv,G*Sigma_r*G',Minv,G*Sigma_r*G',eye(nq)) ...
                        + 4*jacobian_gradient2(mu_r'*G'+Gc',Minv,G*Sigma_r*G',Minv,G*mu_r+Gc) ...
                        + 2*(trace(K*Sigma_r*G'*Jg)+mu_r'*Wx*mu_r+2*mu_r'*Yx+L'*Jg'*Gc) ...
                        *(jacobian_gradient1(Minv,G*Sigma_r*G') + jacobian_gradient(mu_r'*G',Minv,G*mu_r+Gc) ...
                        + jacobian_gradient(mu_r'*G'+Gc',Minv,Gc));
                    
                    dV_M_Drx_Drx_dqdot(i,:) = 0;
                    
                    dV_M_Drx_Dry_dq_first_chain(:,:,i) = -K*Sigma_r*(Wxy+Wxy')*Sigma_r*Ky' - Ky*Sigma_r*Wxy*Sigma_r*K' - K*Sigma_r*Wxy*Sigma_r*Ky' ...
                        -Ky*(mu_r*Oxy'+Oxy*mu_r')*K'-K*(mu_r*Oxy'+Oxy*mu_r')*Ky'-Zy*Oxy'*K'-L*Oxy'*Ky'-K*Oxy*Zy'-Ky*Oxy*L' ...
                        +2*(trace(Ky*Sigma_r*G'*Jg)+mu_r'*Wxy*mu_r+mu_r'*Qxy+Xxy'*mu_r+Zy'*Jg'*Gc) ...
                        *(-K*Sigma_r*Ky' - Ky*mu_r*mu_r'*K' - Ky*mu_r*L' - Zy*mu_r'*K' - Zy*L');
                    dV_M_Drx_Dry_dq(i,:) = trace(dV_M_Drx_Dry_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Drx_Dry_dq(i,:) = dV_M_Drx_Dry_dq(i,:) - 2*E_M_Drx_Dry*dE_M_Drx_Dry_dq(i,:);%[double check this part]
                    dV_M_Drx_Dry_dq(i,:) = dV_M_Drx_Dry_dq(i,:) + jacobian_gradient3(Minv,Fy*Sigma_r*Fy',Minv,G*Sigma_r*G',eye(nq)) ...
                        + jacobian_gradient2(mu_r'*Fy'+Fc',Minv,G*Sigma_r*Fy',Minv,G*mu_r+Gc) ...
                        + jacobian_gradient2(mu_r'*Fy'+Fyc',Minv,G*Sigma_r*G',Minv,Fy*mu_r+Fyc) ...
                        + jacobian_gradient2(mu_r'*G'+Gc',Minv,Fy*Sigma_r*Fy',Minv,G*mu_r+Gc) ...
                        + jacobian_gradient2(mu_r'*G'+Gc',Minv,Fy*Sigma_r*G',Minv,Fy*mu_r+Fyc) ...
                        +2*(trace(Ky*Sigma_r*G'*Jg)+mu_r'*Wxy*mu_r+mu_r'*Qxy+Xxy'*mu_r+Zy'*Jg'*Gc) ...
                        *(jacobian_gradient1(Minv,Fy*Sigma_r*G') + jacobian_gradient(mu_r'*Fy',Minv,G*mu_r) ...
                        + jacobian_gradient(mu_r'*Fy',Minv,Gc) + jacobian_gradient(Fyc',Minv,G*mu_r+Gc));
                    
                    dV_M_Drx_Dry_dqdot(i,:) = 0;
                    
                    dV_M_Dry_Drx_dq(i,:) = dV_M_Drx_Dry_dq(i,:);
                    dV_M_Dry_Drx_dqdot(i,:) = 0;
                end
                
                dE_M_Drx_nr = [0;dE_M_Drx_nr_dq/2;dE_M_Drx_nr_dqdot/2;dE_M_Drx_nr_dq/2;dE_M_Drx_nr_dqdot/2;zeros(3,1);zeros(8,1)];
                dE_M_Drx_Drx = [0;dE_M_Drx_Drx_dq/2;dE_M_Drx_Drx_dqdot/2;dE_M_Drx_Drx_dq/2;dE_M_Drx_Drx_dqdot/2;zeros(3,1);zeros(8,1)];
                dE_M_Drx_Dry = [0;dE_M_Drx_Dry_dq/2;dE_M_Drx_Dry_dqdot/2;dE_M_Drx_Dry_dq/2;dE_M_Drx_Dry_dqdot/2;zeros(3,1);zeros(8,1)];
                dE_M_Dry_Drx = dE_M_Drx_Dry;
                
                % expectation and covariance of b_v_x
                E_b_Drx = (mu_r'*G' + Gc')*J_blk;
                V_b_Drx = trace(T*T'*Sigma_r);
                
                for i=1:nq
                    dE_b_Drx_dq(i,:) = trace( (-h*(K*mu_r + L)*(B*u_prev - C)'*Minv)'*dMdq(:,:,i) ) - trace( h*(K*mu_r + L)'*dCdq(:,i)) ...
                        + trace(((G*mu_r + Gc)*qdot_blk')'*obj.dJgdq_i);% the last part is dJq/dq
                    
                    dE_b_Drx_dqdot(i,:) = - trace( h*(K*mu_r + L)'*dCdqdot(:,i));
                    
                    dV_b_Drx_dq(i,:) = trace( (-h*Minv*(B*u_prev - C)*T'*Sigma_r*K' -h*K*Sigma_r*T*(B*u_prev - C)'*Minv)'*dMdq(:,:,i)) ...
                        - trace( (2*h*K*Sigma_r*T)'*dCdq(:,i)) + jacobian_gradient(G',qdot_blk*qdot_blk',G*Sigma_r);
                    dV_b_Drx_dqdot(i,:) = - trace( (2*h*K*Sigma_r*T)'*dCdqdot(:,i));
                end
                
                dE_b_Drx_dh = (mu_r'*G' + Gc')*Jg*Minv*(B*u_prev - C);
                dE_b_Drx_du = ((mu_r'*G' + Gc')*Jg*Minv*B*h)';
                
                dE_b_Drx = [dE_b_Drx_dh;dE_b_Drx_dq/2;dE_b_Drx_dqdot/2;dE_b_Drx_dq/2;dE_b_Drx_dqdot/2;dE_b_Drx_du;zeros(8,1)];
                
                dV_b_Drx_dh = 2*h*trace(G'*Jg*Minv*(B*u_prev - C)*(B*u_prev - C)'*Minv*Jg'*G*Sigma_r) + trace(G'*Jg*Minv*(B*u_prev - C)*qdot_prev'*Jg'*G*Sigma_r) ...
                    + trace(G'*Jg*qdot_prev*(B*u_prev - C)'*Minv*Jg'*G*Sigma_r);
                dV_b_Drx_du = 2*(G'*Jg*Minv*B)'*h*Sigma_r*T;
                dV_b_Drx = [dV_b_Drx_dh;dV_b_Drx_dq/2;dV_b_Drx_dqdot/2;dV_b_Drx_dq/2;dV_b_Drx_dqdot/2;dV_b_Drx_du;zeros(8,1)];
                
                lambda_n = lambda(1+3*(foot_indx-1));
                lambda_tx = lambda(2+3*(foot_indx-1));
                lambda_ty = lambda(3+3*(foot_indx-1));
                gamma_single = gamma(foot_indx);
                lambda_vec = [lambda_n;lambda_tx;lambda_ty;gamma_single];
                
                E_M_v_x = [h*E_M_Drx_nr, h*E_M_Drx_Drx, h*E_M_Drx_Dry, 1]';
                E_Mvx_lambda_plus_bvx = E_M_v_x'*lambda_vec + E_b_Drx;
                
                dE_M_v_x = [h*dE_M_Drx_nr, h*dE_M_Drx_Drx, h*dE_M_Drx_Dry, zeros(36,1)]';
                
                % NCP residual, currently assume no smoothing func applied
                E_Phi(1+2*(foot_indx-1)) = lambda_tx*E_Mvx_lambda_plus_bvx;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dE_Phi_dh(1+2*(foot_indx-1)) = lambda_tx*(E_M_v_x'*lambda_vec/h + dE_b_Drx_dh);% the last part is dE_b_Drx/dh
                dE_Phi_dq0(:,1+2*(foot_indx-1)) = lambda_tx*h*(dE_M_Drx_nr_dq*lambda_n+dE_M_Drx_Drx_dq*lambda_tx+dE_M_Drx_Dry_dq*lambda_ty)/2 + lambda_tx*dE_b_Drx_dq/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dv0(:,1+2*(foot_indx-1)) = lambda_tx*h*(dE_M_Drx_nr_dqdot*lambda_n+dE_M_Drx_Drx_dqdot*lambda_tx+dE_M_Drx_Dry_dqdot*lambda_ty)/2 + lambda_tx*dE_b_Drx_dqdot/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dq1(:,1+2*(foot_indx-1)) = dE_Phi_dq0(:,1+2*(foot_indx-1));
                dE_Phi_dv1(:,1+2*(foot_indx-1)) = dE_Phi_dv0(:,1+2*(foot_indx-1));
                dE_Phi_du(:,1+2*(foot_indx-1)) = lambda_tx*dE_b_Drx_du;
                dE_Phi_dlambda_n(1+2*(foot_indx-1)) = lambda_tx*h*E_M_Drx_nr;
                dE_Phi_dlambda_tx(1+2*(foot_indx-1)) = E_Mvx_lambda_plus_bvx + h*E_M_Drx_Drx*lambda_tx;
                dE_Phi_dlambda_ty(1+2*(foot_indx-1)) = lambda_tx*h*E_M_Drx_Dry;
                dE_Phi_dgamma(1+2*(foot_indx-1)) = lambda_tx;
                
                if(foot_indx == 1)
                    dE_Phi(:,1) = [dE_Phi_dh(1); dE_Phi_dq0(:,1); dE_Phi_dv0(:,1); dE_Phi_dq1(:,1); dE_Phi_dv1(:,1); dE_Phi_du(:,1);dE_Phi_dlambda_n(1); ...
                        dE_Phi_dlambda_tx(1);dE_Phi_dlambda_ty(1);zeros(3,1);dE_Phi_dgamma(1);0];
                elseif(foot_indx == 2)
                    dE_Phi(:,3) = [dE_Phi_dh(3); dE_Phi_dq0(:,3); dE_Phi_dv0(:,3); dE_Phi_dq1(:,3); dE_Phi_dv1(:,3); dE_Phi_du(:,3);zeros(3,1); ...
                        dE_Phi_dlambda_n(3);dE_Phi_dlambda_tx(3);dE_Phi_dlambda_ty(3);0;dE_Phi_dgamma(3)];
                end
                
                %--------------- third LCP condition ---------------%
                % expectation and covariance of M_v_y
                E_M_Dry_nr = trace(Vy*Sigma_r) + mu_r'*Vy*mu_r + Z'*Jg'*(Fy*mu_r + Fyc) + mu_r'*U'*Jg'*Fyc;
                E_M_Dry_Dry = trace(Wy*Sigma_r) + mu_r'*Wy*mu_r + Ly'*Jg'*(2*Fy*mu_r+Fyc);
                
                for i=1:nq
                    % expectation derivative w.r.t q and qdot
                    dE_M_Dry_nr_dq(i,:) = trace( (-Minv*Jg'*(F*(Sigma_r + mu_r*mu_r')*Fy' + Fc*(Fy*mu_r + Fyc)' + F*mu_r*Fyc')*Jg*Minv)'*dMdq(:,:,i) ) ...
                        + jacobian_gradient(F',Minv,Fy*Sigma_r) + jacobian_gradient(mu_r'*F',Minv,Fy*mu_r) ...
                        + jacobian_gradient(Fc',Minv,Fy*mu_r + Fyc) + jacobian_gradient(mu_r'*F',Minv,Fyc);
                    dE_M_Dry_nr_dqdot(i,:) = 0;
                    
                    dE_M_Dry_Dry_dq(i,:) = trace( (-Minv*Jg'*(Fy*(Sigma_r + mu_r*mu_r')*Fy' + Fyc*(2*Fy*mu_r + Fyc)')*Jg*Minv)'*dMdq(:,:,i) ) ...
                        + jacobian_gradient(Fy',Minv,Fy*Sigma_r) + jacobian_gradient(mu_r'*Fy',Minv,Fy*mu_r) ...
                        + jacobian_gradient(Fyc',Minv,2*Fy*mu_r+Fyc);
                    dE_M_Dry_Dry_dqdot(i,:) = 0;
                    
                    % covariance derivative w.r.t q and qdot
                    dV_M_Dry_nr_dq_first_chain(:,:,i) = -Ky*Sigma_r*(Vy+Vy')*Sigma_r*U' - U*Sigma_r*Vy*Sigma_r*Ky' - Ky*Sigma_r*Vy*Sigma_r*U' ...
                        -U*(mu_r*Oy'+Oy*mu_r')*Ky'-Ky*(mu_r*Oy'+Oy*mu_r')*U'-Z*Oy'*Ky'-Ly*Oy'*U'-Ky*Oy*Z'-U*Oy*Ly' ...
                        +2*(trace(U*Sigma_r*Fy'*Jg)+mu_r'*Vy*mu_r+mu_r'*Qy+Xy'*mu_r+Z'*Jg'*Fyc) ...
                        *(-Ky*Sigma_r*U' - U*mu_r*mu_r'*Ky' - U*mu_r*Ly' - Z*mu_r'*Ky' - Z*Ly');
                    dV_M_Dry_nr_dq_first_chain(:,:,i) = dV_M_Dry_nr_dq_first_chain(:,:,i)';% transponse due to that the equation above is derived based on dV_M_nr_Dry
                    dV_M_Dry_nr_dq(i,:) = trace(dV_M_Dry_nr_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Dry_nr_dq(i,:) = dV_M_Dry_nr_dq(i,:) - 2*E_M_Dry_nr*dE_M_Dry_nr_dq(i,:);%[double check this part]
                    dV_M_Dry_nr_dq(i,:) = dV_M_Dry_nr_dq(i,:) + jacobian_gradient3(Minv,F*Sigma_r*F',Minv,Fy*Sigma_r*Fy',eye(nq)) ...
                        + jacobian_gradient2(mu_r'*F'+Fc',Minv,Fy*Sigma_r*F',Minv,Fy*mu_r+Fyc) ...
                        + jacobian_gradient2(mu_r'*F'+Fc',Minv,Fy*Sigma_r*Fy',Minv,F*mu_r+Fc) ...
                        + jacobian_gradient2(mu_r'*Fy'+Fyc',Minv,F*Sigma_r*F',Minv,Fy*mu_r+Fyc) ...
                        + jacobian_gradient2(mu_r'*Fy'+Fyc',Minv,F*Sigma_r*Fy',Minv,F*mu_r+Fc) ...
                        +2*(trace(U*Sigma_r*Fy'*Jg)+mu_r'*Vy*mu_r+mu_r'*Qy+Xy'*mu_r+Z'*Jg'*Fyc) ...
                        *(jacobian_gradient1(Minv,F*Sigma_r*Fy') + jacobian_gradient(mu_r'*F',Minv,Fy*mu_r) ...
                        + jacobian_gradient(mu_r'*F',Minv,Fyc) + jacobian_gradient(Fc',Minv,Fy*mu_r+Fyc));
                    
                    dV_M_Dry_nr_dqdot(i,:) = 0;
                    
                    dV_M_Dry_Dry_dq_first_chain(:,:,i) = -4*Ky*Sigma_r*Wy*Sigma_r*Ky' + 4*(-Ky*mu_r*mu_r'*Wy*Sigma_r*Ky' - Ky*Sigma_r*Wy*mu_r*mu_r'*Ky' ...
                        -Ky*mu_r*Yy'*Sigma_r*Ky'-Ky*Sigma_r*Wy*mu_r*Ly'-Ly*mu_r'*Wy'*Sigma_r*Ky'-Ky*Sigma_r*Yy*mu_r'*Ky'-Ly*Yy'*Sigma_r*Ky'-Ky*Sigma_r*Yy*Ly') ...
                        +2*(trace(Ky*Sigma_r*Fy'*Jg)+mu_r'*Wy*mu_r+2*mu_r'*Yy+Ly'*Jg'*Fyc) ...
                        *(-Ky*(Sigma_r + mu_r*mu_r')*Ky' - 2*Ky*mu_r*Ly' - Ly*Ly');
                    dV_M_Dry_Dry_dq(i,:) = trace(dV_M_Dry_Dry_dq_first_chain(:,:,i)'*dMdq(:,:,i));
                    dV_M_Dry_Dry_dq(i,:) = dV_M_Dry_Dry_dq(i,:) - 2*E_M_Dry_Dry*dE_M_Dry_Dry_dq(i,:);%[double check this part]
                    dV_M_Dry_Dry_dq(i,:) = dV_M_Dry_Dry_dq(i,:) + 2*jacobian_gradient3(Minv,Fy*Sigma_r*Fy',Minv,Fy*Sigma_r*Fy',eye(nq)) ...
                        + 4*jacobian_gradient2(mu_r'*Fy'+Fyc',Minv,Fy*Sigma_r*Fy',Minv,Fy*mu_r+Fyc) ...
                        + 2*(trace(Ky*Sigma_r*Fy'*Jg)+mu_r'*Wy*mu_r+2*mu_r'*Yy+Ly'*Jg'*Fyc) ...
                        *(jacobian_gradient1(Minv,Fy*Sigma_r*Fy') + jacobian_gradient(mu_r'*Fy',Minv,Fy*mu_r+Fyc) ...
                        + jacobian_gradient(mu_r'*Fy'+Fyc',Minv,Fyc));
                    
                    dV_M_Dry_Dry_dqdot(i,:) = 0;
                end
                
                %[h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dE_M_Dry_nr = [0;dE_M_Dry_nr_dq/2;dE_M_Dry_nr_dqdot/2;dE_M_Dry_nr_dq/2;dE_M_Dry_nr_dqdot/2;zeros(3,1);zeros(8,1)];
                dE_M_Dry_Dry = [0;dE_M_Dry_Dry_dq/2;dE_M_Dry_Dry_dqdot/2;dE_M_Dry_Dry_dq/2;dE_M_Dry_Dry_dqdot/2;zeros(3,1);zeros(8,1)];
                
                % expectation and covariance of b_v_y
                E_b_Dry = (mu_r'*Fy' + Fyc')*J_blk;
                V_b_Dry = trace(Ty*Ty'*Sigma_r);
                
                for i=1:nq
                    dE_b_Dry_dq(i,:) = trace( (-h*(Ky*mu_r + Ly)*(B*u_prev - C)'*Minv)'*dMdq(:,:,i) ) - trace( h*(Ky*mu_r + Ly)'*dCdq(:,i)) ...
                        + trace(((Fy*mu_r + Fyc)*qdot_blk')'*obj.dJgdq_i);% the last part is dJq/dq;
                    dE_b_Dry_dqdot(i,:) = - trace( h*(Ky*mu_r + Ly)'*dCdqdot(:,i));
                    
                    dV_b_Dry_dq(i,:) = trace( (-h*Minv*(B*u_prev - C)*Ty'*Sigma_r*Ky' -h*Ky*Sigma_r*Ty*(B*u_prev - C)'*Minv)'*dMdq(:,:,i)) ...
                        - trace( (2*h*Ky*Sigma_r*Ty)'*dCdq(:,i)) + jacobian_gradient(Fy',qdot_blk*qdot_blk',Fy*Sigma_r);
                    dV_b_Dry_dqdot(i,:) = - trace( (2*h*Ky*Sigma_r*Ty)'*dCdqdot(:,i));
                end
                dE_b_Dry_dh = (mu_r'*Fy' + Fyc')*Jg*Minv*(B*u_prev - C);
                dE_b_Dry_du = ((mu_r'*Fy' + Fyc')*Jg*Minv*B*h)';
                
                dE_b_Dry = [dE_b_Dry_dh;dE_b_Dry_dq/2;dE_b_Dry_dqdot/2;dE_b_Dry_dq/2;dE_b_Dry_dqdot/2;dE_b_Dry_du;zeros(8,1)];
                
                dV_b_Dry_dh = 2*h*trace(Fy'*Jg*Minv*(B*u_prev - C)*(B*u_prev - C)'*Minv*Jg'*Fy*Sigma_r) + trace(Fy'*Jg*Minv*(B*u_prev - C)*qdot_prev'*Jg'*Fy*Sigma_r) ...
                    + trace(Fy'*Jg*qdot_prev*(B*u_prev - C)'*Minv*Jg'*Fy*Sigma_r);
                dV_b_Dry_du = 2*(Fy'*Jg*Minv*B)'*h*Sigma_r*Ty;
                dV_b_Dry = [dV_b_Dry_dh;dV_b_Dry_dq/2;dV_b_Dry_dqdot/2;dV_b_Dry_dq/2;dV_b_Dry_dqdot/2;dV_b_Dry_du;zeros(8,1)];
                
                % vectrozie expectation components
                E_M_v_y = [h*E_M_Dry_nr, h*E_M_Dry_Drx, h*E_M_Dry_Dry, 1]';
                E_Mvy_lambda_plus_bvy = E_M_v_y'*lambda_vec + E_b_Dry;
                
                dE_M_v_y = [h*dE_M_Dry_nr, h*dE_M_Dry_Drx, h*dE_M_Dry_Dry, zeros(36,1)]';
                
                % NCP residual, currently assume no smoothing func applied
                E_Phi(2+2*(foot_indx-1)) = lambda_ty*E_Mvy_lambda_plus_bvy;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dE_Phi_dh(2+2*(foot_indx-1)) = lambda_ty*(E_M_v_y'*lambda_vec/h + dE_b_Dry_dh);% the last part is dE_b_Dry/dh
                dE_Phi_dq0(:,2+2*(foot_indx-1)) = lambda_ty*h*(dE_M_Dry_nr_dq*lambda_n+dE_M_Dry_Drx_dq*lambda_tx+dE_M_Dry_Dry_dq*lambda_ty)/2 + lambda_ty*dE_b_Dry_dq/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dv0(:,2+2*(foot_indx-1)) = lambda_ty*h*(dE_M_Dry_nr_dqdot*lambda_n+dE_M_Dry_Drx_dqdot*lambda_tx+dE_M_Dry_Dry_dqdot*lambda_ty)/2 + lambda_ty*dE_b_Dry_dqdot/2;% the last 1/2 is due to d((q0+q1)/2)/dq0
                dE_Phi_dq1(:,2+2*(foot_indx-1)) = dE_Phi_dq0(:,2+2*(foot_indx-1));
                dE_Phi_dv1(:,2+2*(foot_indx-1)) = dE_Phi_dv0(:,2+2*(foot_indx-1));
                dE_Phi_du(:,2+2*(foot_indx-1)) = lambda_ty*dE_b_Dry_du;
                dE_Phi_dlambda_n(2+2*(foot_indx-1)) = lambda_ty*h*E_M_Dry_nr;
                dE_Phi_dlambda_tx(2+2*(foot_indx-1)) = lambda_ty*h*E_M_Dry_Drx;
                dE_Phi_dlambda_ty(2+2*(foot_indx-1)) = E_Mvy_lambda_plus_bvy + lambda_ty*h*E_M_Dry_Dry;
                dE_Phi_dgamma(2+2*(foot_indx-1)) = lambda_ty;
                
                if(foot_indx == 1)
                    dE_Phi(:,2) = [dE_Phi_dh(2); dE_Phi_dq0(:,2); dE_Phi_dv0(:,2); dE_Phi_dq1(:,2); dE_Phi_dv1(:,2); dE_Phi_du(:,2);dE_Phi_dlambda_n(2); ...
                        dE_Phi_dlambda_tx(2);dE_Phi_dlambda_ty(2);zeros(3,1);dE_Phi_dgamma(2);0];
                elseif(foot_indx == 2)
                    dE_Phi(:,4) = [dE_Phi_dh(4); dE_Phi_dq0(:,4); dE_Phi_dv0(:,4); dE_Phi_dq1(:,4); dE_Phi_dv1(:,4); dE_Phi_du(:,4);zeros(3,1); ...
                        dE_Phi_dlambda_n(4);dE_Phi_dlambda_tx(4);dE_Phi_dlambda_ty(4);0;dE_Phi_dgamma(4)];
                end
                
                % LCP variance matrix of V_Mvx_lambda_plus_bvx
                %fourth order expectation
                [E_xnxn,dE_xnxn_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,F,Fc);
                [E_xnxx,dE_xnxx_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,G,Gc);
                [E_xnxy,dE_xnxy_dq] = expectation_fourth_order_multiply(G,Gc,F,Fc,G,Gc,H,Hc);
                
                [E_xxxn,dE_xxxn_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,F,Fc);
                [E_xxxx,dE_xxxx_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,G,Gc);
                [E_xxxy,dE_xxxy_dq] = expectation_fourth_order_multiply(G,Gc,G,Gc,G,Gc,H,Hc);
                
                [E_xyxn,dE_xyxn_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,F,Fc);
                [E_xyxx,dE_xyxx_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,G,Gc);
                [E_xyxy,dE_xyxy_dq] = expectation_fourth_order_multiply(G,Gc,H,Hc,G,Gc,H,Hc);
                
                E_Mvx_lambda_lambda_Mvx_quad = h^2*(E_xnxn*lambda_n^2 + E_xnxx*lambda_n*lambda_tx + E_xnxy*lambda_n*lambda_ty ...
                    + E_xxxn*lambda_tx*lambda_n + E_xxxx*lambda_tx^2 + E_xxxy*lambda_tx*lambda_ty ...
                    + E_xyxn*lambda_ty*lambda_n + E_xyxx*lambda_ty*lambda_tx + E_xyxy*lambda_ty^2);
                E_Mvx_lambda_lambda_Mvx_linear = 2*h*E_M_Drx_nr*lambda_n*gamma_single + 2*h*E_M_Drx_Drx*lambda_tx*gamma_single ...
                    + 2*h*E_M_Drx_Dry*lambda_ty*gamma_single + gamma_single^2;
                E_Mvx_lambda_lambda_Mvx = E_Mvx_lambda_lambda_Mvx_quad + E_Mvx_lambda_lambda_Mvx_linear;
                
                % computing derivative of E_Mvx_lambda_lambda_Mvx
                dE_Mvx_lambda_lambda_Mvx_dh = 2*E_Mvx_lambda_lambda_Mvx_quad/h + (E_Mvx_lambda_lambda_Mvx_linear-gamma_single^2)/h;
                dE_Mvx_lambda_lambda_Mvx_dq0 = h^2/2*(dE_xnxn_dq*lambda_n^2 + dE_xnxx_dq*lambda_n*lambda_tx + dE_xnxy_dq*lambda_n*lambda_ty ...
                    + dE_xxxn_dq*lambda_tx*lambda_n + dE_xxxx_dq*lambda_tx^2 + dE_xxxy_dq*lambda_tx*lambda_ty ...
                    + dE_xyxn_dq*lambda_ty*lambda_n + dE_xyxx_dq*lambda_ty*lambda_tx + dE_xyxy_dq*lambda_ty^2) + h*dE_M_Drx_nr_dq*lambda_n*gamma_single ...
                    + h*dE_M_Drx_Drx_dq*lambda_tx*gamma_single + h*dE_M_Drx_Dry_dq*lambda_ty*gamma_single;
                dE_Mvx_lambda_lambda_Mvx_dv0 = zeros(nv,1);
                dE_Mvx_lambda_lambda_Mvx_dq1 = dE_Mvx_lambda_lambda_Mvx_dq0;
                dE_Mvx_lambda_lambda_Mvx_dv1 = zeros(nv,1);
                dE_Mvx_lambda_lambda_Mvx_du = zeros(nu,1);
                dE_Mvx_lambda_lambda_Mvx_dlambda_n = h^2*(2*E_xnxn*lambda_n + E_xnxx*lambda_tx + E_xnxy*lambda_ty + E_xxxn*lambda_tx + E_xyxn*lambda_ty) ...
                    + 2*h*E_M_Drx_nr*gamma_single;
                dE_Mvx_lambda_lambda_Mvx_dlambda_tx = h^2*(E_xnxx*lambda_n + E_xxxn*lambda_n + 2*E_xxxx*lambda_tx + E_xxxy*lambda_ty + E_xyxx*lambda_ty) ...
                    + 2*h*E_M_Drx_Drx*gamma_single;
                dE_Mvx_lambda_lambda_Mvx_dlambda_ty = h^2*(E_xnxy*lambda_n + E_xxxy*lambda_tx + E_xyxn*lambda_n + E_xyxx*lambda_tx + 2* E_xyxy*lambda_ty) ...
                    + 2*h*E_M_Drx_Dry*gamma_single;
                dE_Mvx_lambda_lambda_Mvx_dgamma = 2*h*E_M_Drx_nr*lambda_n + 2*h*E_M_Drx_Drx*lambda_tx ...
                    + 2*h*E_M_Drx_Dry*lambda_ty + 2*gamma_single;
                
                if(foot_indx == 1)
                    dE_Mvx_lambda_lambda_Mvx = [dE_Mvx_lambda_lambda_Mvx_dh;dE_Mvx_lambda_lambda_Mvx_dq0;dE_Mvx_lambda_lambda_Mvx_dv0;dE_Mvx_lambda_lambda_Mvx_dq1; ...
                        dE_Mvx_lambda_lambda_Mvx_dv1;dE_Mvx_lambda_lambda_Mvx_du;dE_Mvx_lambda_lambda_Mvx_dlambda_n;dE_Mvx_lambda_lambda_Mvx_dlambda_tx; ...
                        dE_Mvx_lambda_lambda_Mvx_dlambda_ty;zeros(3,1);dE_Mvx_lambda_lambda_Mvx_dgamma;0];
                elseif(foot_indx == 2)
                    dE_Mvx_lambda_lambda_Mvx = [dE_Mvx_lambda_lambda_Mvx_dh;dE_Mvx_lambda_lambda_Mvx_dq0;dE_Mvx_lambda_lambda_Mvx_dv0;dE_Mvx_lambda_lambda_Mvx_dq1; ...
                        dE_Mvx_lambda_lambda_Mvx_dv1;dE_Mvx_lambda_lambda_Mvx_du;zeros(3,1);dE_Mvx_lambda_lambda_Mvx_dlambda_n; ...
                        dE_Mvx_lambda_lambda_Mvx_dlambda_tx;dE_Mvx_lambda_lambda_Mvx_dlambda_ty;0;dE_Mvx_lambda_lambda_Mvx_dgamma];
                end
                
                %cov(M_{v,x}^T,b_{v,x}) = E[(M_{v,x}^T - E(M_{v,x}^T))(b_{v,x}-E(b_{v,x}))] = E[M_{v,x}^T*b_{v,x}] - E[M_{v,x}^T]*E[b_{v,x}];
                % E[M_{v,x}^T*b_{v,x}]
                [E_Mvx_bvx_Drx_nr_Drx,dE_Mvx_bvx_Drx_nr_Drx_dq,dE_Mvx_bvx_Drx_nr_Drx_dqdot,dE_Mvx_bvx_Drx_nr_Drx_dh,dE_Mvx_bvx_Drx_nr_Drx_du] = expectation_third_order_multiply(G,Gc,F,Fc,G,Gc);
                [E_Mvx_bvx_Drx_Drx_Drx,dE_Mvx_bvx_Drx_Drx_Drx_dq,dE_Mvx_bvx_Drx_Drx_Drx_dqdot,dE_Mvx_bvx_Drx_Drx_Drx_dh,dE_Mvx_bvx_Drx_Drx_Drx_du] = expectation_third_order_multiply(G,Gc,G,Gc,G,Gc);
                [E_Mvx_bvx_Drx_Dry_Drx,dE_Mvx_bvx_Drx_Dry_Drx_dq,dE_Mvx_bvx_Drx_Dry_Drx_dqdot,dE_Mvx_bvx_Drx_Dry_Drx_dh,dE_Mvx_bvx_Drx_Dry_Drx_du] = expectation_third_order_multiply(G,Gc,Fy,Fyc,G,Gc);
                E_Mvx_bvx = [E_Mvx_bvx_Drx_nr_Drx;E_Mvx_bvx_Drx_Drx_Drx;E_Mvx_bvx_Drx_Dry_Drx;E_b_Drx];
                
                dE_Mvx_bvx_Drx_nr_Drx = [dE_Mvx_bvx_Drx_nr_Drx_dh;dE_Mvx_bvx_Drx_nr_Drx_dq/2;dE_Mvx_bvx_Drx_nr_Drx_dqdot/2;dE_Mvx_bvx_Drx_nr_Drx_dq/2;dE_Mvx_bvx_Drx_nr_Drx_dqdot/2;dE_Mvx_bvx_Drx_nr_Drx_du;zeros(8,1)];
                dE_Mvx_bvx_Drx_Drx_Drx = [dE_Mvx_bvx_Drx_Drx_Drx_dh;dE_Mvx_bvx_Drx_Drx_Drx_dq/2;dE_Mvx_bvx_Drx_Drx_Drx_dqdot/2;dE_Mvx_bvx_Drx_Drx_Drx_dq/2;dE_Mvx_bvx_Drx_Drx_Drx_dqdot/2;dE_Mvx_bvx_Drx_Drx_Drx_du;zeros(8,1)];
                dE_Mvx_bvx_Drx_Dry_Drx = [dE_Mvx_bvx_Drx_Dry_Drx_dh;dE_Mvx_bvx_Drx_Dry_Drx_dq/2;dE_Mvx_bvx_Drx_Dry_Drx_dqdot/2;dE_Mvx_bvx_Drx_Dry_Drx_dq/2;dE_Mvx_bvx_Drx_Dry_Drx_dqdot/2;dE_Mvx_bvx_Drx_Dry_Drx_du;zeros(8,1)];
                
                dE_Mvx_bvx = [dE_Mvx_bvx_Drx_nr_Drx,dE_Mvx_bvx_Drx_Drx_Drx,dE_Mvx_bvx_Drx_Dry_Drx,dE_b_Drx]';
                
                %  V[M_{v,x}*lambda+b_{v,x}] = V[M_{v,x}*lambda] + cov(M_{v,x}*lambda,b_{v,x}) + cov(b_{v,x},M_{v,x}*lambda) + V[b_{v,x}]
                % = V[M_{v,x}*lambda] + 2*lambda^T*cov(M_{v,x}^T,b_{v,x}) + V[b_{v,x}]
                % = E[M_{v,x}*lambda*lambda^T*M_{v,x}^T] - (E[M_{v,x}*lambda])^2 + 2*lambda^T*cov(M_{v,x}^T,b_{v,x}) + V[b_{v,x}]
                % = E[M_{v,x}*lambda*lambda^T*M_{v,x}^T] - (E[M_{v,x}*lambda])^2 + 2*lambda^T*(E[M_{v,x}^T*b_{v,x}] - E[M_{v,x}^T]*E[b_{v,x}]) + V[b_{v,x}]
                V_Mvx_lambda_plus_bvx = E_Mvx_lambda_lambda_Mvx - (E_M_v_x'*lambda_vec)^2 + 2*lambda_vec'*(E_Mvx_bvx - E_M_v_x*E_b_Drx) + V_b_Drx;%[double check]
                
                V_Phi(1+2*(foot_indx-1)) = lambda_tx^2*V_Mvx_lambda_plus_bvx;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dV_Phi(:,1+2*(foot_indx-1)) = lambda_tx^2*(dE_Mvx_lambda_lambda_Mvx - 2*(E_M_v_x'*lambda_vec)*(dE_M_v_x'*lambda_vec) + (2*lambda_vec'*dE_Mvx_bvx)' - 2*lambda_vec'*E_M_v_x*dE_b_Drx ...
                    -(2*lambda_vec'*dE_M_v_x*E_b_Drx)' + dV_b_Drx);
                
                % LCP variance matrix of V_Mvy_lambda_plus_bvy
                %fourth order expectation
                [E_ynyn,dE_ynyn_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,F,Fc);
                [E_ynyx,dE_ynyx_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,G,Gc);
                [E_ynyy,dE_ynyy_dq] = expectation_fourth_order_multiply(H,Hc,F,Fc,H,Hc,H,Hc);
                
                [E_yxyn,dE_yxyn_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,F,Fc);
                [E_yxyx,dE_yxyx_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,G,Gc);
                [E_yxyy,dE_yxyy_dq] = expectation_fourth_order_multiply(H,Hc,G,Gc,H,Hc,H,Hc);
                
                [E_yyyn,dE_yyyn_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,F,Fc);
                [E_yyyx,dE_yyyx_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,G,Gc);
                [E_yyyy,dE_yyyy_dq] = expectation_fourth_order_multiply(H,Hc,H,Hc,H,Hc,H,Hc);
                
                % computing E_Mvy_lambda_lambda_Mvy
                E_Mvy_lambda_lambda_Mvy_quad = h^2*(E_ynyn*lambda_n^2 + E_ynyx*lambda_n*lambda_tx + E_ynyy*lambda_n*lambda_ty ...
                    + E_yxyn*lambda_tx*lambda_n + E_yxyx*lambda_tx^2 + E_yxyy*lambda_tx*lambda_ty ...
                    + E_yyyn*lambda_ty*lambda_n + E_yyyx*lambda_ty*lambda_tx + E_yyyy*lambda_ty^2);
                E_Mvy_lambda_lambda_Mvy_linear = 2*h*E_M_Dry_nr*lambda_n*gamma_single + 2*h*E_M_Dry_Drx*lambda_tx*gamma_single ...
                    + 2*h*E_M_Dry_Dry*lambda_ty*gamma_single + gamma_single^2;
                E_Mvy_lambda_lambda_Mvy = E_Mvy_lambda_lambda_Mvy_quad + E_Mvy_lambda_lambda_Mvy_linear;
                
                % computing derivative of E_Mvy_lambda_lambda_Mvy
                dE_Mvy_lambda_lambda_Mvy_dh = 2*E_Mvy_lambda_lambda_Mvy_quad/h + (E_Mvy_lambda_lambda_Mvy_linear-gamma_single^2)/h;
                dE_Mvy_lambda_lambda_Mvy_dq0 = h^2/2*(dE_ynyn_dq*lambda_n^2 + dE_ynyx_dq*lambda_n*lambda_tx + dE_ynyy_dq*lambda_n*lambda_ty ...
                    + dE_yxyn_dq*lambda_tx*lambda_n + dE_yxyx_dq*lambda_tx^2 + dE_yxyy_dq*lambda_tx*lambda_ty ...
                    + dE_yyyn_dq*lambda_ty*lambda_n + dE_yyyx_dq*lambda_ty*lambda_tx + dE_yyyy_dq*lambda_ty^2) + h*dE_M_Dry_nr_dq*lambda_n*gamma_single ...
                    + h*dE_M_Dry_Drx_dq*lambda_tx*gamma_single + h*dE_M_Dry_Dry_dq*lambda_ty*gamma_single;
                dE_Mvy_lambda_lambda_Mvy_dv0 = zeros(nv,1);
                dE_Mvy_lambda_lambda_Mvy_dq1 = dE_Mvy_lambda_lambda_Mvy_dq0;
                dE_Mvy_lambda_lambda_Mvy_dv1 = zeros(nv,1);
                dE_Mvy_lambda_lambda_Mvy_du = zeros(nu,1);
                dE_Mvy_lambda_lambda_Mvy_dlambda_n = h^2*(2*E_ynyn*lambda_n + E_ynyx*lambda_tx + E_ynyy*lambda_ty + E_yxyn*lambda_tx + E_yyyn*lambda_ty) ...
                    + 2*h*E_M_Dry_nr*gamma_single;
                dE_Mvy_lambda_lambda_Mvy_dlambda_tx = h^2*(E_ynyx*lambda_n + E_yxyn*lambda_n + 2*E_yxyx*lambda_tx + E_yxyy*lambda_ty + E_yyyx*lambda_ty) ...
                    + 2*h*E_M_Dry_Drx*gamma_single;
                dE_Mvy_lambda_lambda_Mvy_dlambda_ty = h^2*(E_ynyy*lambda_n + E_yxyy*lambda_tx + E_yyyn*lambda_n + E_yyyx*lambda_tx + 2* E_yyyy*lambda_ty) ...
                    + 2*h*E_M_Dry_Dry*gamma_single;
                dE_Mvy_lambda_lambda_Mvy_dgamma = 2*h*E_M_Dry_nr*lambda_n + 2*h*E_M_Dry_Drx*lambda_tx ...
                    + 2*h*E_M_Dry_Dry*lambda_ty + 2*gamma_single;
                
                if(foot_indx == 1)
                    dE_Mvy_lambda_lambda_Mvy = [dE_Mvy_lambda_lambda_Mvy_dh;dE_Mvy_lambda_lambda_Mvy_dq0;dE_Mvy_lambda_lambda_Mvy_dv0;dE_Mvy_lambda_lambda_Mvy_dq1; ...
                        dE_Mvy_lambda_lambda_Mvy_dv1;dE_Mvy_lambda_lambda_Mvy_du;dE_Mvy_lambda_lambda_Mvy_dlambda_n; ...
                        dE_Mvy_lambda_lambda_Mvy_dlambda_tx;dE_Mvy_lambda_lambda_Mvy_dlambda_ty;zeros(3,1);dE_Mvy_lambda_lambda_Mvy_dgamma;0];
                elseif(foot_indx == 2)
                    dE_Mvy_lambda_lambda_Mvy = [dE_Mvy_lambda_lambda_Mvy_dh;dE_Mvy_lambda_lambda_Mvy_dq0;dE_Mvy_lambda_lambda_Mvy_dv0;dE_Mvy_lambda_lambda_Mvy_dq1; ...
                        dE_Mvy_lambda_lambda_Mvy_dv1;dE_Mvy_lambda_lambda_Mvy_du;zeros(3,1);dE_Mvy_lambda_lambda_Mvy_dlambda_n; ...
                        dE_Mvy_lambda_lambda_Mvy_dlambda_tx;dE_Mvy_lambda_lambda_Mvy_dlambda_ty;0;dE_Mvy_lambda_lambda_Mvy_dgamma];
                end
                %Computing cov(M_{v,y}^T,b_{v,y}) = E[(M_{v,y}^T - E(M_{v,y}^T))(b_{v,y}-E(b_{v,y}))] = E[M_{v,y}^T*b_{v,y}] - E[M_{v,y}^T]*E[b_{v,y}];
                
                % E[M_{v,y}^T*b_{v,y}]
                [E_Mvy_bvy_Dry_nr_Dry,dE_Mvy_bvy_Dry_nr_Dry_dq,dE_Mvy_bvy_Dry_nr_Dry_dqdot,dE_Mvy_bvy_Dry_nr_Dry_dh,dE_Mvy_bvy_Dry_nr_Dry_du] = expectation_third_order_multiply(H,Hc,F,Fc,H,Hc);
                [E_Mvy_bvy_Dry_Drx_Dry,dE_Mvy_bvy_Dry_Drx_Dry_dq,dE_Mvy_bvy_Dry_Drx_Dry_dqdot,dE_Mvy_bvy_Dry_Drx_Dry_dh,dE_Mvy_bvy_Dry_Drx_Dry_du] = expectation_third_order_multiply(H,Hc,G,Gc,H,Hc);
                [E_Mvy_bvy_Dry_Dry_Dry,dE_Mvy_bvy_Dry_Dry_Dry_dq,dE_Mvy_bvy_Dry_Dry_Dry_dqdot,dE_Mvy_bvy_Dry_Dry_Dry_dh,dE_Mvy_bvy_Dry_Dry_Dry_du] = expectation_third_order_multiply(H,Hc,H,Hc,H,Hc);
                E_Mvy_bvy = [E_Mvy_bvy_Dry_nr_Dry;E_Mvy_bvy_Dry_Drx_Dry;E_Mvy_bvy_Dry_Dry_Dry;E_b_Dry];
                
                dE_Mvy_bvy_Dry_nr_Dry = [dE_Mvy_bvy_Dry_nr_Dry_dh;dE_Mvy_bvy_Dry_nr_Dry_dq/2;dE_Mvy_bvy_Dry_nr_Dry_dqdot/2;dE_Mvy_bvy_Dry_nr_Dry_dq/2;dE_Mvy_bvy_Dry_nr_Dry_dqdot/2;dE_Mvy_bvy_Dry_nr_Dry_du;zeros(8,1)];
                dE_Mvy_bvy_Dry_Drx_Dry = [dE_Mvy_bvy_Dry_Drx_Dry_dh;dE_Mvy_bvy_Dry_Drx_Dry_dq/2;dE_Mvy_bvy_Dry_Drx_Dry_dqdot/2;dE_Mvy_bvy_Dry_Drx_Dry_dq/2;dE_Mvy_bvy_Dry_Drx_Dry_dqdot/2;dE_Mvy_bvy_Dry_Drx_Dry_du;zeros(8,1)];
                dE_Mvy_bvy_Dry_Dry_Dry = [dE_Mvy_bvy_Dry_Dry_Dry_dh;dE_Mvy_bvy_Dry_Dry_Dry_dq/2;dE_Mvy_bvy_Dry_Dry_Dry_dqdot/2;dE_Mvy_bvy_Dry_Dry_Dry_dq/2;dE_Mvy_bvy_Dry_Dry_Dry_dqdot/2;dE_Mvy_bvy_Dry_Dry_Dry_du;zeros(8,1)];
                
                dE_Mvy_bvy = [dE_Mvy_bvy_Dry_nr_Dry,dE_Mvy_bvy_Dry_Drx_Dry,dE_Mvy_bvy_Dry_Dry_Dry,dE_b_Dry]';
                
                %  V[M_{v,y}*lambda+b_{v,y}] = V[M_{v,y}*lambda] + cov(M_{v,y}*lambda,b_{v,y}) + cov(b_{v,y},M_{v,y}*lambda) + V[b_{v,y}]
                % = V[M_{v,y}*lambda] + 2*lambda^T*cov(M_{v,y}^T,b_{v,y}) + V[b_{v,y}]
                % = E[M_{v,y}*lambda*lambda^T*M_{v,y}^T] - (E[M_{v,y}*lambda])^2 + 2*lambda^T*cov(M_{v,y}^T,b_{v,y}) + V[b_{v,y}]
                % = E[M_{v,y}*lambda*lambda^T*M_{v,y}^T] - (E[M_{v,y}*lambda])^2 + 2*lambda^T*(E[M_{v,y}^T*b_{v,y}] - E[M_{v,y}^T]*E[b_{v,y}]) + V[b_{v,y}]
                V_Mvy_lambda_plus_bvy = E_Mvy_lambda_lambda_Mvy - (E_M_v_y'*lambda_vec)^2 + 2*lambda_vec'*(E_Mvy_bvy - E_M_v_y*E_b_Dry) + V_b_Dry;
                
                V_Phi(2+2*(foot_indx-1)) = lambda_ty^2*V_Mvy_lambda_plus_bvy;
                % derivative w.r.t [h;q0;q0dot;q1;q1dot;u;lambda_n;lambda_tx;lambda_ty;gamma]
                dV_Phi(:,2+2*(foot_indx-1)) = lambda_ty^2*(dE_Mvy_lambda_lambda_Mvy - 2*(E_M_v_y'*lambda_vec)*(dE_M_v_y'*lambda_vec) + (2*lambda_vec'*dE_Mvy_bvy)' - 2*lambda_vec'*E_M_v_y*dE_b_Dry ...
                    -(2*lambda_vec'*dE_M_v_y*E_b_Dry)'+ dV_b_Dry);
            end
            
            delta = 10^7;% coefficient
            f = delta/2 * (norm(E_Phi)^2 + norm(V_Phi)^2);% - slack_var*ones(zdim,1);
            df = delta*(E_Phi'*dE_Phi' + V_Phi'*dV_Phi');
            %df = delta*(gamma.*g)'*[diag(gamma)*dg + [zeros(zdim,xdim) diag(g)]];
            
            persistent sliding_LCP_ERM_NCP_residual
            sliding_LCP_ERM_NCP_residual = [sliding_LCP_ERM_NCP_residual, f/(delta/2)];
            if length(sliding_LCP_ERM_NCP_residual) == obj.N-1
                sliding_LCP_ERM_NCP_residual
                sliding_LCP_ERM_NCP_residual = [];
            end
            
            % obj.verbose_print = 1;
            % if obj.verbose_print == 1
            %     disp('ERM NCP residual square');
            %     f
            % end
            
            %% debugging
            %------- begin comparison ----------
            % a comparison test between probabilistic expectation and distribution-free Phi
            v1_est = v0 + Minv*(B*u_prev - C + D{1}'*[lambda(2);lambda(5)] + D{2}'*[lambda(3);lambda(6)] + n'*[lambda(1);lambda(4)])*h;
            
            %deterministic (distribution-free) cost func and its derivative
            f_deter = zeros(obj.nC*(1+obj.nD),1);
            df_deter = zeros(obj.nC*(1+obj.nD),nq+nv+obj.nC*(2+obj.nD));
            
            f_deter(1:1+obj.nD:end) = phi;
            df_deter(1:1+obj.nD:end,1:nq) = n;
            for j=1:obj.nD
                f_deter(1+j:1+obj.nD:end) = gamma+D{j}*v1;
                df_deter(1+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                df_deter(1+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                df_deter(1+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v1);%d/dq
            end
            
            Phi = f_deter.*lambda;
            %E_Phi
            %Phi
            %% ------- end comparison ----------
            
            % % test non-zero values
            % nonzero_index_set = find(abs(dE_Phi) > 1e-3);
            % if length(nonzero_index_set) > 4
            %     disp('number of nonzero index set elements > 4')
            % end
            
            function dX_dq = jacobian_gradient(A,B,C)
                %X = trace(A*Jg*B*Jg'*C)
                dX_dJg = A'*C'*obj.Jg*B' + C*A*obj.Jg*B;
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient1(A,B)
                %X = trace(A*Jg'*B*Jg)
                dX_dJg = B*obj.Jg*A + B'*obj.Jg*A';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient2(A,B,C,D,E)
                %X = trace(A*Jg*B*Jg'*C*Jg*D*Jg'*E)
                dX_dJg = A'*E'*obj.Jg*D'*obj.Jg'*C'*obj.Jg*B' + C*obj.Jg*D*obj.Jg'*E*A*obj.Jg*B ...
                    + C'*obj.Jg*B'*obj.Jg'*A'*E'*obj.Jg*D' + E*A*obj.Jg*B*obj.Jg'*C*obj.Jg*D;
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient3(A,B,C,D,E)
                %X = trace(A*Jg'*B*Jg*C*Jg'*D*Jg*E)
                dX_dJg = B*obj.Jg*C*obj.Jg'*D*obj.Jg*E*A + B'*obj.Jg*A'*E'*obj.Jg'*D'*obj.Jg*C' ...
                    + D*obj.Jg*E*A*obj.Jg'*B*obj.Jg*C + D'*obj.Jg*C'*obj.Jg'*B'*obj.Jg*A'*E';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient4(A,B,C,D,E)
                %X = trace(A*Jg'*B*Jg)*C*Jg*D*Jg'*E
                dX_dJg = B*obj.Jg*A*(C*obj.Jg*D*obj.Jg'*E) + B'*obj.Jg*A'*(E'*obj.Jg*D'*obj.Jg'*C') ...
                    + trace(A*Jg'*B*Jg)*C'*E'*obj.Jg*D' + trace(A*Jg'*B*Jg)*E*C*obj.Jg*D;
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient5(A,B,C,D,E,F)
                %X = A*Jg*B*Jg'*C*trace(D*Jg'*E*Jg*F)
                dX_dJg = A'*C'*obj.Jg*B'*trace(D*Jg'*E*Jg*F) + C*A*obj.Jg*B*trace(D*Jg'*E*Jg*F) ...
                    + A*Jg*B*Jg'*C*E*obj.Jg*F*D + A*Jg*B*Jg'*C*E'*obj.Jg*D'*F';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient6(A,B,C,D)
                %X = A*Jg*B*Jg'*C*Jg*D
                dX_dJg = A'*D'*obj.Jg'*C'*obj.Jg*B' + C*obj.Jg*D*A*obj.Jg*B + C'*obj.Jg*B'*obj.Jg'*A'*D';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient7(A,B,C,D)
                %X = trace(A*Jg'*B*Jg)*C*Jg*D
                dX_dJg = B*obj.Jg*A*(C*obj.Jg*D) + B'*obj.Jg*A'*(C*obj.Jg*D) + trace(A*obj.Jg'*B*obj.Jg)*C'*D';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function dX_dq = jacobian_gradient8(A,B,C,D)
                %X = A*Jg*B*Jg'*C*Jg*D
                dX_dJg = A'*D'*obj.Jg'*C'*obj.Jg*B' + C*obj.Jg*D*A*obj.Jg*B + C'*obj.Jg*B'*obj.Jg'*A'*D';
                
                dX_dq = trace(dX_dJg'*obj.dJgdq_i);
            end
            
            function [E,dEdq,dEdqdot,dEdh,dEdu] = expectation_third_order_multiply(Ain, ain, Bin, bin, Cin, cin)
                %refactor input matrices
                AAin = obj.h*obj.Minv*obj.Jg'*Ain;
                aain = obj.h*obj.Minv*obj.Jg'*ain;
                BBin = obj.Jg'*Bin;
                bbin = obj.Jg'*bin;
                CCin = obj.J_blk'*Cin;
                ccin = obj.J_blk'*cin;
                
                % expectation
                term1 = ((AAin*obj.mu_r+aain)'*BBin+(BBin*obj.mu_r+bbin)'*AAin);
                term2 = CCin';
                term3 = trace(AAin*obj.Sigma_r*BBin')+(AAin*obj.mu_r+aain)'*(BBin*obj.mu_r+bbin);
                term4 = (CCin*obj.mu_r+ccin)';
                
                E = term1*obj.Sigma_r*term2 +term3*term4;
                
                % derivative of expectation (d/dM * dM/dq)
                % first term
                dEdM = - (AAin*obj.mu_r+aain)*term2'*obj.Sigma_r*BBin'*obj.Minv - obj.Minv*(BBin*obj.mu_r+bbin)*term2'*obj.Sigma_r*AAin' - obj.Minv*obj.h*obj.Jg'*Cin*obj.Sigma_r*term1'*(obj.B*obj.u-C)'*obj.Minv;
                % second term
                dEdM = dEdM -obj.Minv*term4*BBin*obj.Sigma_r*AAin' -(AAin*obj.mu_r+aain)*term4*(BBin*obj.mu_r+bbin)'*obj.Minv - obj.Minv*(obj.Jg'*Cin*obj.mu_r+obj.Jg'*cin)*term3*obj.h*(obj.B*obj.u-obj.C)'*obj.Minv;
                
                dEdC = -(term1*Sigma_r*Cin'*Jg*h*Minv)' - (term3*(mu_r'*Cin'+cin')*Jg*Minv*h)';
                
                dEdh = E/obj.h + term1*obj.Sigma_r*Cin'*obj.Jg*obj.Minv*(obj.B*obj.u - obj.C) +term3*(obj.mu_r'*Cin' + cin')*obj.Jg*obj.Minv*(obj.B*obj.u - obj.C);
                
                dEdJ_times_dJdq = jacobian_gradient6(obj.mu_r'*Ain'+ain',obj.Minv*obj.h,Bin*obj.Sigma_r*Cin',obj.qdot_blk) ...
                    + jacobian_gradient6(obj.mu_r'*Bin'+bin',obj.h*obj.Minv,Ain*obj.Sigma_r*Cin',obj.qdot_blk) ...
                    + jacobian_gradient7(obj.h*obj.Minv,Ain*obj.Sigma_r*Bin',obj.mu_r'*Cin',obj.qdot_blk) ...
                    + jacobian_gradient8(obj.mu_r'*Ain'+ain',obj.Minv*obj.h,(Bin*obj.mu_r+bin)*obj.mu_r'*Cin',obj.qdot_blk);
                
                for i=1:nq
                    dEdq(i,:) = trace(dEdM'*obj.dMdq(:,:,i)) + trace(dEdC'*obj.dCdq(:,i)) + dEdJ_times_dJdq;
                    dEdqdot(i,:) = trace(dEdC'*obj.dCdqdot(:,i));
                end
                
                dEdu = (term1*obj.Sigma_r*Cin'*obj.Jg*(obj.Minv*obj.B*obj.h))' - (term3*(obj.mu_r'*Cin'+cin')*obj.Jg*obj.Minv*obj.B*obj.h)';
            end
            
            function [E,dEdq] = expectation_fourth_order_multiply(Ain, ain, Bin, bin, Cin, cin, Din, din)
                %refactor input matrices
                AAin = obj.Minv*obj.Jg'*Ain;
                aain = obj.Minv*obj.Jg'*ain;
                BBin = obj.Jg'*Bin;
                bbin = obj.Jg'*bin;
                CCin = obj.Minv*obj.Jg'*Cin;
                ccin = obj.Minv*obj.Jg'*cin;
                DDin = obj.Jg'*Din;
                ddin = obj.Jg'*din;
                
                % expectation
                term1 = ((AAin*obj.mu_r+aain)'*BBin+(BBin*obj.mu_r+bbin)'*AAin);
                term2 = (CCin'*(DDin*obj.mu_r+ddin)+DDin'*(CCin*obj.mu_r+ccin));
                term3 = trace(AAin*obj.Sigma_r*BBin')+(AAin*obj.mu_r+aain)'*(BBin*obj.mu_r+bbin);
                term4 = trace(CCin*obj.Sigma_r*DDin')+(CCin*obj.mu_r+ccin)'*(DDin*obj.mu_r+ddin);
                
                E = trace(AAin*obj.Sigma_r*(CCin'*DDin + DDin'*CCin)*obj.Sigma_r*BBin') + term1*obj.Sigma_r*term2 + term3*term4;
                
                % derivative of expectation (d/dM * dM/dq)
                % first term
                dEdM = -obj.Minv*BBin*obj.Sigma_r*(CCin'*DDin+DDin'*CCin)*obj.Sigma_r*AAin' - CCin*obj.Sigma_r*AAin'*BBin*obj.Sigma_r*DDin'*obj.Minv - obj.Minv*DDin*obj.Sigma_r*AAin'*BBin*obj.Sigma_r*CCin';
                % second term
                dEdM = dEdM - obj.Minv*(AAin*obj.mu_r+aain)*term2'*obj.Sigma_r*BBin'*obj.Minv - obj.Minv*(BBin*obj.mu_r+bbin)*term2'*obj.Sigma_r*AAin' - CCin*obj.Sigma_r*term1'*(DDin*obj.mu_r+ddin)'*obj.Minv ...
                    - obj.Minv*DDin*obj.Sigma_r*term1'*(obj.mu_r'*CCin'+ccin');
                % third term
                dEdM = dEdM -obj.Minv*term4*BBin*obj.Sigma_r*AAin' -(AAin*obj.mu_r+aain)*term4*(BBin*obj.mu_r+bbin)'*obj.Minv - obj.Minv*term3*DDin*obj.Sigma_r*CCin'...
                    - (CCin*obj.mu_r+ccin)*term3*(DDin*obj.mu_r+ddin)'*obj.Minv;
                
                dEdJ_times_dJdq = jacobian_gradient3(obj.Minv,Ain*obj.Sigma_r*Cin',obj.Minv,Din*obj.Sigma_r*Bin',eye(nq)) ...
                    + jacobian_gradient3(obj.Minv,Ain*obj.Sigma_r*Din',obj.Minv,Cin*obj.Sigma_r*Bin',eye(nq)) ...
                    + jacobian_gradient2(obj.mu_r'*Ain'+ain',obj.Minv,Bin*obj.Sigma_r*Cin',obj.Minv,Din*obj.mu_r+din) ...
                    + jacobian_gradient2(obj.mu_r'*Ain'+ain',obj.Minv,Bin*obj.Sigma_r*Din',obj.Minv,Cin*obj.mu_r+cin) ...
                    + jacobian_gradient2(obj.mu_r'*Bin'+bin',obj.Minv,Ain*obj.Sigma_r*Cin',obj.Minv,Din*obj.mu_r+din) ...
                    + jacobian_gradient2(obj.mu_r'*Bin'+bin',obj.Minv,Ain*obj.Sigma_r*Din',obj.Minv,Cin*obj.mu_r+cin) ...
                    + jacobian_gradient3(obj.Minv,Ain*obj.Sigma_r*Bin',obj.Minv,Cin*obj.Sigma_r*Din',eye(obj.nq)) ...
                    + jacobian_gradient2(obj.mu_r'*Ain'+ain',obj.Minv,(Bin*obj.mu_r+bin)*(obj.mu_r'*Cin'+cin'),obj.Minv,Din*obj.mu_r+din) ...
                    + jacobian_gradient4(obj.Minv,Ain*obj.Sigma_r*Bin',obj.mu_r'*Cin'+cin',obj.Minv,Din*obj.mu_r+din) ...
                    + jacobian_gradient5(obj.mu_r'*Ain'+ain',obj.Minv,(Bin*obj.mu_r+bin),obj.Minv,Cin*obj.Sigma_r*Din',eye(obj.nq));
                
                %dEdC = 0;
                
                for i=1:nq
                    dEdq(i,:) = trace(dEdM'*obj.dMdq(:,:,i)) + dEdJ_times_dJdq;
                    dEdqot(i,:) = 0;
                end
                
            end
            
        end
        
        function [f,df] = dynamics_constraint_fun(obj,h,x0,x1,u,lambda,lambda_jl,Fext)
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nu = obj.plant.getNumInputs;
            nl = length(lambda);
            njl = length(lambda_jl);
            obj.h = h;
            
            lambda = lambda*obj.options.lambda_mult;
            lambda_jl = lambda_jl*obj.options.lambda_jl_mult;
            
            assert(nq == nv) % not quite ready for the alternative
            
            q0 = x0(1:nq);
            v0 = x0(nq+1:nq+nv);
            q1 = x1(1:nq);
            v1 = x1(nq+1:nq+nv);
            
            switch obj.options.integration_method
                case RobustContactImplicitTrajectoryOptimization_Brick.MIDPOINT
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics((q0+q1)/2,(v0+v1)/2);
                    dH0 = dH/2;
                    dC0 = dC/2;
                    dB0 = dB/2;
                    dH1 = dH/2;
                    dC1 = dC/2;
                    dB1 = dB/2;
                case RobustContactImplicitTrajectoryOptimization_Brick.FORWARD_EULER
                    [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                    dH1 = zeros(nq^2,2*nq);
                    dC1 = zeros(nq,2*nq);
                    dB1 = zeros(nq*nu,2*nq);
                case RobustContactImplicitTrajectoryOptimization_Brick.BACKWARD_EULER
                    [H,C,B,dH1,dC1,dB1] = obj.plant.manipulatorDynamics(q1,v1);
                    dH0 = zeros(nq^2,2*nq);
                    dC0 = zeros(nq,2*nq);
                    dB0 = zeros(nq*nu,2*nq);
                case RobustContactImplicitTrajectoryOptimization_Brick.MIXED
                    [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                    dH1 = zeros(nq^2,2*nq);
                    dC1 = zeros(nq,2*nq);
                    dB1 = zeros(nq*nu,2*nq);
            end
            
            BuminusC = B*u-C + [Fext;zeros(5,1)];% manually add an external force
            if nu>0,
                dBuminusC0 = matGradMult(dB0,u) - dC0;
                dBuminusC1 = matGradMult(dB1,u) - dC1;
            else
                dBuminusC0 = -dC0;
                dBuminusC1 = -dC1;
            end
            
            switch obj.options.integration_method
                case RobustContactImplicitTrajectoryOptimization_Brick.MIDPOINT
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*(v0 + v1)/2;
                    dfq = [-(v1+v0)/2, -eye(nq), -h/2*eye(nq), eye(nq), -h/2*eye(nq) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization_Brick.FORWARD_EULER
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*v0;
                    dfq = [-v0, -eye(nq), -h*eye(nq), eye(nq), zeros(nq,nv) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization_Brick.BACKWARD_EULER
                    % q1 = q0 + h*v1
                    fq = q1 - q0 - h*v1;
                    dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq) zeros(nq,nu+nl+njl)];
                case RobustContactImplicitTrajectoryOptimization_Brick.MIXED
                    fq = q1 - q0 - h*v1;
                    %dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl)];
                    dfqdFext = zeros(nq,1);
                    dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl), dfqdFext];% expand the dimension with Fext
            end
            
            % H*v1 = H*v0 + h*(B*u - C) + n^T lambda_N + d^T * lambda_f
            fv = H*(v1 - v0) - h*BuminusC;
            % [h q0 v0 q1 v1 u l ljl]
            
            %dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl)] + ...
            %    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,nu+nl+njl)];
            % expand the dimension with Fext
            dfvdFext = -h*[1;zeros(5,1)];
            dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl), dfvdFext] + ...
                [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,1+nu+nl+njl)];
            
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
            
            f_numeric = f;
            df_numeric = df;
            %            disp('check gradient')
            
            %             [f_numeric,df_numeric] = geval(@(h,x0,x1,lambda,Fext) dynamics_constraint_fun_check(obj,h,x0,x1,u,lambda,lambda_jl,Fext),h,x0,x1,lambda,Fext,struct('grad_method','numerical'));
            %             valuecheck(df,df_numeric,1e-3);
            %             valuecheck(f,f_numeric,1e-3);
            
            function [f_num,df_num] = dynamics_constraint_fun_check(obj,h,x0,x1,u,lambda,lambda_jl,Fext)
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
                    case RobustContactImplicitTrajectoryOptimization_Brick.MIDPOINT
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics((q0+q1)/2,(v0+v1)/2);
                        dH0 = dH/2;
                        dC0 = dC/2;
                        dB0 = dB/2;
                        dH1 = dH/2;
                        dC1 = dC/2;
                        dB1 = dB/2;
                    case RobustContactImplicitTrajectoryOptimization_Brick.FORWARD_EULER
                        [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                        dH1 = zeros(nq^2,2*nq);
                        dC1 = zeros(nq,2*nq);
                        dB1 = zeros(nq*nu,2*nq);
                    case RobustContactImplicitTrajectoryOptimization_Brick.BACKWARD_EULER
                        [H,C,B,dH1,dC1,dB1] = obj.plant.manipulatorDynamics(q1,v1);
                        dH0 = zeros(nq^2,2*nq);
                        dC0 = zeros(nq,2*nq);
                        dB0 = zeros(nq*nu,2*nq);
                    case RobustContactImplicitTrajectoryOptimization_Brick.MIXED
                        [H,C,B,dH0,dC0,dB0] = obj.plant.manipulatorDynamics(q0,v0);
                        dH1 = zeros(nq^2,2*nq);
                        dC1 = zeros(nq,2*nq);
                        dB1 = zeros(nq*nu,2*nq);
                end
                
                BuminusC = B*u-C + [Fext;zeros(5,1)];% manually add an external force
                if nu>0,
                    dBuminusC0 = matGradMult(dB0,u) - dC0;
                    dBuminusC1 = matGradMult(dB1,u) - dC1;
                else
                    dBuminusC0 = -dC0;
                    dBuminusC1 = -dC1;
                end
                
                switch obj.options.integration_method
                    case RobustContactImplicitTrajectoryOptimization_Brick.MIDPOINT
                        % q1 = q0 + h*v1
                        fq = q1 - q0 - h*(v0 + v1)/2;
                        dfq = [-(v1+v0)/2, -eye(nq), -h/2*eye(nq), eye(nq), -h/2*eye(nq) zeros(nq,nu+nl+njl)];
                    case RobustContactImplicitTrajectoryOptimization_Brick.FORWARD_EULER
                        % q1 = q0 + h*v1
                        fq = q1 - q0 - h*v0;
                        dfq = [-v0, -eye(nq), -h*eye(nq), eye(nq), zeros(nq,nv) zeros(nq,nu+nl+njl)];
                    case RobustContactImplicitTrajectoryOptimization_Brick.BACKWARD_EULER
                        % q1 = q0 + h*v1
                        fq = q1 - q0 - h*v1;
                        dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq) zeros(nq,nu+nl+njl)];
                    case RobustContactImplicitTrajectoryOptimization_Brick.MIXED
                        fq = q1 - q0 - h*v1;
                        %dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl)];
                        dfqdFext = zeros(nq,1);
                        dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl), dfqdFext];% expand the dimension with Fext
                end
                
                % H*v1 = H*v0 + h*(B*u - C) + n^T lambda_N + d^T * lambda_f
                fv = H*(v1 - v0) - h*BuminusC;
                % [h q0 v0 q1 v1 u l ljl]
                
                %dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl)] + ...
                %    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,nu+nl+njl)];
                % expand the dimension with Fext
                dfvdFext = [1;zeros(5,1)];
                dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl), dfvdFext] + ...
                    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,1+nu+nl+njl)];
                
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
                
                f_num = [fq;fv];
                df_num = [dfq;dfv];
            end
            
        end
        
        function [xtraj,utraj,ltraj,ljltraj,slacktraj,Fexttraj,z,F,info,infeasible_constraint_name] = solveTraj(obj,t_init,traj_init)
            [xtraj,utraj,z,F,info,infeasible_constraint_name] = solveTraj@DirectTrajectoryOptimization_Brick(obj,t_init,traj_init);
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
            Fexttraj = PPTrajectory(foh(t,[reshape(z(obj.F_ext_inds),1, obj.N)]));
        end
        
        function obj = setupVariables(obj,N)
            obj = setupVariables@DirectTrajectoryOptimization_Brick(obj,N);
            [~,normal,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), false, obj.options.active_collision_options);
            obj.nC = size(normal, 2);
            obj.nD = 2*length(d);
            obj.nFext = 1;
            
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
            
            obj.F_ext_inds = reshape(obj.num_vars + (1:N), 1, N);
            obj = obj.addDecisionVariable(N);
            
            obj.LCP_slack_inds = reshape(obj.num_vars + (1:N-1), 1, N-1);
            obj = obj.addDecisionVariable(N-1);
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
            end
            
            % initialize external force variables
            if isfield(traj_init,'F_ext')
                F_ext_eval = traj_init.F_ext.eval(t_init);
                z0(obj.F_ext_inds) = F_ext_eval(:,1:obj.N); % take N values
            else
                z0(obj.F_ext_inds) = 0;
            end
            
            if obj.nC > 0
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
            nFext = obj.nFext;
            running_cost = FunctionHandleObjective(1+nX+nFext,running_cost_function);
            for i=1:obj.N-1
                obj = obj.addCost(running_cost,{obj.h_inds(i);obj.x_inds(:,i);obj.F_ext_inds(:,i)});
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
                    
                    obj.nonlincompl_constraints_purturb{i} = NonlinearComplementarityConstraint_original(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode,obj.options.compl_slack);
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
                    
                    lincompl_constraints{i} = LinearComplementarityConstraint_original(W,r,M,obj.options.lincc_mode,obj.options.lincompl_slack);
                    obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds]);
                end
                
                if obj.nJL > 0
                    % joint limit linear complementarity constraint
                    % lambda_jl /perp [q - lb_jl; -q + ub_jl]
                    W_jl = zeros(obj.nJL);
                    [r_jl,M_jl] = jointLimitConstraints(plant_perturb,q0);
                    jlcompl_constraints{i} = LinearComplementarityConstraint_original(W_jl,r_jl,M_jl,obj.options.lincc_mode,obj.options.jlcompl_slack);
                    
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