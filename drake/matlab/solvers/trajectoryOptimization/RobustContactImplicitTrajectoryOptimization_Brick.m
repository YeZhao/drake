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
        MIXED = 4;   % matched to TimeSteppingRigidBodyManipulator_Brick. Forward on qd, backward on q
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
            global uncertainty_source;
            
            constraints = cell(N-1,1);
            lincompl_constraints = cell(N-1,1);
            obj.nonlincompl_constraints = cell(N-1,1);
            obj.nonlincompl_slack_inds = cell(N-1,1);
            jlcompl_constraints = cell(N-1,1);
            dyn_inds = cell(N-1,1);
            
            n_vars = 2*nX + nU + 1 + obj.nC*(2+obj.nD) + obj.nJL + obj.nFext;
            cnstr = FunctionHandleConstraint(zeros(nX,1),zeros(nX,1),n_vars,@obj.dynamics_constraint_fun);
            q0 = getZeroConfiguration(obj.plant);
            
            %[~,~,~,~,~,~,~,mu] = obj.plant.contactConstraints(q0,false,obj.options.active_collision_options);
            manip = obj.plant.getManipulator();
            if strcmp(manip.uncertainty_source, 'friction_coeff') || strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                mu = mean(manip.uncertain_mu_set)*ones(obj.nC,1);
            else
                mu = ones(obj.nC,1);%nominal value
            end
            
            if strcmp(manip.uncertainty_source, 'friction_coeff')
                uncertainty_source = 'friction_coeff';
            elseif strcmp(manip.uncertainty_source, 'terrain_height')
                uncertainty_source = 'terrain_height';
            elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                uncertainty_source = 'friction_coeff+terrain_height';
            elseif isempty(manip.uncertainty_source)
                uncertainty_source = [];
            end
            
            external_force_max = 5;% random chosen value
            external_force_cnstr = BoundingBoxConstraint([-external_force_max,-external_force_max],[external_force_max,external_force_max]);
            
            %contact_force_cnstr = BoundingBoxConstraint(zeros(obj.nC*(1+obj.nD),1),inf(obj.nC*(1+obj.nD),1));
            contact_force_cnstr = BoundingBoxConstraint(-1*ones(obj.nC*(1+obj.nD),1),inf(obj.nC*(1+obj.nD),1));
            sliding_velocity_cnstr = BoundingBoxConstraint(-1*ones(obj.nC,1),inf(obj.nC,1));
            
            for i=1:obj.N-1
                dyn_inds{i} = {obj.h_inds(i);obj.x_inds(:,i);obj.x_inds(:,i+1);obj.u_inds(:,i);obj.l_inds(:,i);obj.ljl_inds(:,i);obj.F_ext_inds(:,i)};
                constraints{i} = cnstr.setName(sprintf('dynamics[%d]',i));
                %obj = obj.addConstraint(constraints{i}, dyn_inds{i});
                
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
                    
                    %obj.options.nlcc_mode = 5;% robust mode
                    %obj.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode);
                    %obj.nonlincompl_slack_inds{i} = obj.num_vars+1:obj.num_vars + obj.nonlincompl_constraints{i}.n_slack;
                    %obj = obj.addConstraint(obj.nonlincompl_constraints{i},[obj.x_inds(:,i+1);gamma_inds;lambda_inds;obj.LCP_slack_inds(:,i)]);
                    
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
                    
                    %lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode);
                    %obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds;obj.LCP_slack_inds(:,i)]);
                    
                    % ERM formulation
                    assert(size(lambda_inds,1) == obj.nC*(1+obj.nD));
                    assert(size(gamma_inds,1) == obj.nC);
                    
                    % more direct way to implement lambda and gamma nonnegativity
                    % contact_force_bounding_constraints{i} = contact_force_cnstr.setName(sprintf('contactForceConstraint[%d]',i));
                    % obj = obj.addConstraint(contact_force_bounding_constraints{i}, [lambda_inds]);
                    %
                    % sliding_velocity_bounding_constraints{i} = sliding_velocity_cnstr.setName(sprintf('slidingVelocityConstraint[%d]',i));
                    % obj = obj.addConstraint(sliding_velocity_bounding_constraints{i}, [gamma_inds]);
                    
                    % % nonlinear LCP non-negative constraint
                    %obj.options.nlcc_mode = 7;% lambda non-negative mode
                    %obj.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + obj.nC,obj.nC*(1+obj.nD),obj.options.nlcc_mode);
                    %obj = obj.addConstraint(obj.nonlincompl_constraints{i},[obj.x_inds(:,i+1);gamma_inds;lambda_inds;obj.LCP_slack_inds(:,i)]);
                    %
                    % % linear LCP non-negative constraint
                    %obj.options.lincc_mode = 4;% tangential velocity non-negative mode
                    %lincompl_constraints{i} = LinearComplementarityConstraint(W,r,M,obj.options.lincc_mode);
                    %obj = obj.addConstraint(lincompl_constraints{i},[lambda_inds;gamma_inds;obj.LCP_slack_inds(:,i)]);
                    %
                    %inner product residual minimization
                    %obj = obj.addCost(FunctionHandleObjective(nX + obj.nC+obj.nC*(1+obj.nD),@(x,gamma,lambda)ERMcost(obj,x,gamma,lambda),1), ...
                    %    {obj.x_inds(:,i+1);gamma_inds;lambda_inds});
                    
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
            
            % a hacky way to obtain updated timestep (no cost is added, just want to get time step h)
            obj = obj.addCost(FunctionHandleObjective(1,@(h_inds)getTimeStep(obj,h_inds),1),{obj.h_inds(1)});
             
            % robust variance cost with state feedback control
            x_inds_stack = reshape(obj.x_inds,obj.N*nX,[]);
            u_inds_stack = reshape(obj.u_inds,obj.N*nU,[]);
            Fext_inds_stack = reshape(obj.F_ext_inds,obj.N*nFext,[]);
            lambda_inds_stack = reshape(obj.lambda_inds,(obj.N-1)*nL,[]);
            obj.cached_Px = zeros(obj.nx,obj.nx,obj.N);
            obj.cached_Px(:,:,1) = obj.options.Px_coeff*eye(obj.nx); %[ToDo: To be modified]
            %obj = obj.addCost(FunctionHandleObjective(obj.N*(nX+nU),@(x_inds,Fext_inds)robustVariancecost_ML(obj,x_inds,Fext_inds),1),{x_inds_stack;Fext_inds_stack});
            obj = obj.addCost(FunctionHandleObjective(obj.N*(nX+nU),@(x_inds,Fext_inds)robustVariancecost_scaled(obj,x_inds,Fext_inds),1),{x_inds_stack;Fext_inds_stack});
            
            %if (obj.nC > 0)
            %    obj = obj.addCost(FunctionHandleObjective(length(obj.LCP_slack_inds),@(slack)robustLCPcost(obj,slack),1),obj.LCP_slack_inds(:));
            %end
            
            global uncertain_mu;
            global terrain_index;
            
            function [f,df] = ERMcost(obj,x,gamma,lambda)
                y = [x;gamma;lambda];
                global f_accumulate;
                global number_of_call;
                
                if strcmp(manip.uncertainty_source, 'friction_coeff')
                    sample_num = length(manip.uncertain_mu_set);
                    for i=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(i);
                        uncertain_mu = manip.uncertain_mu_set(i);
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif strcmp(manip.uncertainty_source, 'terrain_height')
                    sample_num = length(manip.uncertain_position_set);
                    for i=1:sample_num
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        terrain_index = i;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                    sample_num = length(manip.uncertain_mu_set);
                    for i=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        uncertain_mu = manip.uncertain_mu_set(i);
                        terrain_index = i;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif isempty(manip.uncertainty_source)%non-robust version
                    sample_num = 1;
                    [ff,dff] = ERM_nonlincompl_fun(obj,y);
                    ff_stack = ff;
                    dff_stack = dff;
                else
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,58))*sqrt(0.01);
                    %y = (1-2*rand(1,58))*sqrt(0.01);
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                end
                
                ff_mean = sum(ff_stack,2)/sample_num;
                ff_mean_vec = repmat(ff_mean,1,sample_num);
                dff_mean = sum(dff_stack,3)/sample_num;
                mean_dev = zeros(1,length(y));
                var_dev = zeros(1,length(y));
                
                decision_vec = [gamma;lambda];
                
                for i=1:sample_num
                    mean_dev = mean_dev + decision_vec'*ff_stack(:,i)*decision_vec'*dff_stack(:,:,i);
                    var_dev = var_dev + (decision_vec'*(ff_stack(:,i) - ff_mean))'*(decision_vec'*dff_stack(:,:,i) - decision_vec'*dff_mean);
                end
                
                delta = obj.options.ERMcost_coeff;% coefficient
                f = delta/2 * (norm(decision_vec'*ff_stack)^2 + norm(decision_vec'*ff_stack - decision_vec'*ff_mean_vec)^2);
                df = delta*(mean_dev + var_dev) ...
                    + [zeros(1,nX), delta*((decision_vec'*ff_stack)*ff_stack' + (decision_vec'*(ff_stack - ff_mean_vec))*(ff_stack - ff_mean_vec)')];
                
                if isempty(f_accumulate)
                    f_accumulate = f;
                    number_of_call = 1;
                elseif (number_of_call == obj.N-2)
                    fprintf('ERM cost: %4.4f\n',f);
                    f_accumulate = [];
                    number_of_call = 0;
                else
                    f_accumulate = f_accumulate + f;
                    number_of_call = number_of_call + 1;
                end
                
                f_numeric = f;
                df_numeric = df;
                X0 = [x;gamma;lambda];
                
                % [f_numeric,df_numeric] = geval(@(X0) ERMcost_check(X0),X0,struct('grad_method','numerical'));
                % valuecheck(df,df_numeric,1e-5);
                % valuecheck(f,f_numeric,1e-5);
                
                function [f,df] = ERMcost_check(X0)
                    x = X0(1:obj.nx);
                    gamma = X0(obj.nx+1:obj.nx+obj.nC);
                    lambda = X0(obj.nx+obj.nC+1:end);
                    y = X0;
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        sample_num = length(manip.plant.uncertain_mu_set);
                        for i=1:sample_num
                            manip.uncertain_mu = manip.uncertain_mu_set(i);
                            uncertain_mu = manip.uncertain_mu_set(i);
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        sample_num = length(manip.uncertain_position_set);
                        for i=1:sample_num
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            terrain_index = i;
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        sample_num = length(manip.uncertain_mu_set);
                        for i=1:sample_num
                            manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            uncertain_mu = manip.uncertain_mu_set(i);
                            terrain_index = i;
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif isempty(manip.uncertainty_source)%non-robust version
                        sample_num = 1;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack = ff;
                        dff_stack = dff;
                    else
                        w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                        save -ascii friction_coeff_noise.dat w_mu
                    end
                    
                    ff_mean = sum(ff_stack,2)/sample_num;
                    ff_mean_vec = repmat(ff_mean,1,sample_num);
                    dff_mean = sum(dff_stack,3)/sample_num;
                    mean_dev = zeros(1,length(y));
                    var_dev = zeros(1,length(y));
                    
                    decision_vec = [gamma;lambda];
                    
                    for i=1:sample_num
                        mean_dev = mean_dev + decision_vec'*ff_stack(:,i)*decision_vec'*dff_stack(:,:,i);
                        var_dev = var_dev + (decision_vec'*(ff_stack(:,i) - ff_mean))'*(decision_vec'*dff_stack(:,:,i) - decision_vec'*dff_mean);
                    end
                    
                    delta = obj.options.ERMcost_coeff;% coefficient
                    f = delta/2 * (norm(decision_vec'*ff_stack)^2 + norm(decision_vec'*ff_stack - decision_vec'*ff_mean_vec)^2);
                    df = delta*(mean_dev + var_dev) ...
                        + [zeros(1,nX), delta*((decision_vec'*ff_stack)*ff_stack' + (decision_vec'*(ff_stack - ff_mean_vec))*(ff_stack - ff_mean_vec)')];
                end
            end
            
            function [f,df] = ERMcost_NCP(obj,x,gamma,lambda)
                y = [x;gamma;lambda];
                global f_accumulate;
                global number_of_call;
                global NCP_slack_param
                
                if strcmp(manip.uncertainty_source, 'friction_coeff')
                    sample_num = length(manip.uncertain_mu_set);
                    for i=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(i);
                        uncertain_mu = manip.uncertain_mu_set(i);
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif strcmp(manip.uncertainty_source, 'terrain_height')
                    sample_num = length(manip.uncertain_position_set);
                    for i=1:sample_num
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        terrain_index = i;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                    sample_num = length(manip.uncertain_mu_set);
                    for i=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        uncertain_mu = manip.uncertain_mu_set(i);
                        terrain_index = i;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack(:,i) = ff;
                        dff_stack(:,:,i) = dff;
                    end
                elseif isempty(manip.uncertainty_source)%non-robust version
                    sample_num = 1;
                    [ff,dff] = ERM_nonlincompl_fun(obj,y);
                    ff_stack = ff;
                    dff_stack = dff;
                else
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,58))*sqrt(0.01);
                    %y = (1-2*rand(1,58))*sqrt(0.01);
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                end
                
                ff_mean = sum(ff_stack,2)/sample_num;
                ff_mean_vec = repmat(ff_mean,1,sample_num);
                dff_mean = sum(dff_stack,3)/sample_num;
                mean_dev = zeros(1,length(y));
                var_dev = zeros(1,length(y));
                
                decision_vec = [gamma;lambda];
                decision_length = length(decision_vec);
                
                delta = obj.options.ERMcost_coeff;% coefficient
                
                f = 0;
                df = zeros(1,length(y));
                s = NCP_slack_param;
                for i=1:sample_num
                    mean_residual = decision_vec + ff_stack(:,i) - sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s);
                    var_residual = decision_vec + (ff_stack(:,i)-ff_mean) - sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s);
                    f = f + delta/2 * (norm(mean_residual)^2 + norm(var_residual)^2);
                    mean_dev = mean_residual'*(dff_stack(:,:,i) - diag(ff_stack(:,i)./sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s))*dff_stack(:,:,i));
                    var_dev = var_residual'*((dff_stack(:,:,i) - dff_mean) - diag((ff_stack(:,i)-ff_mean)./sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s))*(dff_stack(:,:,i) - dff_mean));
                    df = df + delta*(mean_dev+var_dev);
                    df = df + [zeros(1,nX), delta*(mean_residual'*(eye(decision_length) - diag(decision_vec./sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s))))];
                    df = df + [zeros(1,nX), delta*(var_residual'*(eye(decision_length) - diag(decision_vec./sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s))))];
                end
                
                if isempty(f_accumulate)
                    f_accumulate = f;
                    number_of_call = 1;
                elseif (number_of_call == obj.N-2)
                    fprintf('ERM cost: %4.4f\n',f);
                    f_accumulate = [];
                    number_of_call = 0;
                else
                    f_accumulate = f_accumulate + f;
                    number_of_call = number_of_call + 1;
                end
                
                %f_numeric = f;
                %df_numeric = df;
                %X0 = [x;gamma;lambda];
                
                % [f_numeric,df_numeric] = geval(@(X0) ERMcost_NCP_check(X0),X0,struct('grad_method','numerical'));
                % valuecheck(df,df_numeric,1e-5);
                % valuecheck(f,f_numeric,1e-5);
                
                function [f,df] = ERMcost_NCP_check(X0)
                    x = X0(1:obj.nx);
                    gamma = X0(obj.nx+1:obj.nx+obj.nC);
                    lambda = X0(obj.nx+obj.nC+1:end);
                    y = X0;
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        sample_num = length(manip.uncertain_mu_set);
                        for i=1:sample_num
                            manip.uncertain_mu = manip.uncertain_mu_set(i);
                            uncertain_mu = manip.uncertain_mu_set(i);
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        sample_num = length(manip.uncertain_position_set);
                        for i=1:sample_num
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            terrain_index = i;
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        sample_num = length(manip.uncertain_mu_set);
                        for i=1:sample_num
                            manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            uncertain_mu = manip.uncertain_mu_set(i);
                            terrain_index = i;
                            [ff,dff] = ERM_nonlincompl_fun(obj,y);
                            ff_stack(:,i) = ff;
                            dff_stack(:,:,i) = dff;
                        end
                    elseif isempty(manip.uncertainty_source)%non-robust version
                        sample_num = 1;
                        [ff,dff] = ERM_nonlincompl_fun(obj,y);
                        ff_stack = ff;
                        dff_stack = dff;
                    else
                        w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                        save -ascii friction_coeff_noise.dat w_mu
                        %x = (1-2*rand(1,58))*sqrt(0.01);
                        %y = (1-2*rand(1,58))*sqrt(0.01);
                        %w_phi = [x;y];%height noise
                        %save -ascii initial_position_noise.dat w_phi
                    end
                    
                    ff_mean = sum(ff_stack,2)/sample_num;
                    ff_mean_vec = repmat(ff_mean,1,sample_num);
                    dff_mean = sum(dff_stack,3)/sample_num;
                    mean_dev = zeros(1,length(y));
                    var_dev = zeros(1,length(y));
                    
                    decision_vec = [gamma;lambda];
                    decision_length = length(decision_vec);
                    
                    delta = obj.options.ERMcost_coeff;% coefficient
                    
                    f = 0;
                    df = zeros(1,length(y));
                    s = 1e-10;
                    for i=1:sample_num
                        mean_residual = decision_vec + ff_stack(:,i) - sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s);
                        var_residual = decision_vec + (ff_stack(:,i)-ff_mean) - sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s);
                        f = f + delta/2 * (norm(mean_residual)^2 + norm(var_residual)^2);
                        mean_dev = mean_residual'*(dff_stack(:,:,i) - diag(ff_stack(:,i)./sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s))*dff_stack(:,:,i));
                        var_dev = var_residual'*((dff_stack(:,:,i) - dff_mean) - diag((ff_stack(:,i)-ff_mean)./sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s))*(dff_stack(:,:,i) - dff_mean));
                        df = df + delta*(mean_dev+var_dev);
                        df = df + [zeros(1,nX), delta*(mean_residual'*(eye(decision_length) - diag(decision_vec./sqrt(decision_vec.^2 + ff_stack(:,i).^2 + s))))];
                        df = df + [zeros(1,nX), delta*(var_residual'*(eye(decision_length) - diag(decision_vec./sqrt(decision_vec.^2 + (ff_stack(:,i)-ff_mean).^2 + s))))];
                    end
                end
            end
            
            % nonlinear complementarity constraints for ERM formulation
            % the only difference here is that the first input is obj since we
            % we have to use obj membership parameters for uncertainty expression:
            %   lambda_N /perp phi(q)
            %   lambda_fi /perp gamma + Di*psi(q,v)
            % x = [q;v;gamma]
            % z = [lambda_N;lambda_F1;lambda_f2] (each contact sequentially)
            function [f,df] = ERM_nonlincompl_fun(obj,y)
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                x = y(1:nq+nv+obj.nC);
                z = y(nq+nv+obj.nC+1:end);
                gamma = x(nq+nv+1:end);
                q = x(1:nq);
                v = x(nq+1:nq+nv);
                
                if strcmp(manip.uncertainty_source, 'friction_coeff')
                    [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                elseif strcmp(manip.uncertainty_source, 'terrain_height')
                    [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                    if isempty(terrain_index)
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    else
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                    end
                elseif isempty(manip.uncertainty_source)
                    [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                end
                
                %% debugging
                % a test of using q0 to derive v1, instead of using v1
                % x0 = zeros(12,1);
                % q0 = x0(1:nq);%zeros(12,1);
                % v0 = x0(nq+1:nq+nv);
                % h = 0.01;
                % u = zeros(3,1);
                % [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                %
                % [M,C,B,dM,dC,dB] = obj.plant.manipulatorDynamics(q,v);
                % Minv = inv(M);
                %
                %v1_est = v0 + Minv*(B*u - C + D{1}'*[z(2);z(5)] + D{2}'*[z(3);z(6)] + n'*[z(1);z(4)])*h;
                %% end debugging
                
                f = zeros(obj.nC*(2+obj.nD),1);
                df = zeros(obj.nC*(2+obj.nD),nq+nv+obj.nC*(2+obj.nD));
                
                %add linear complementarity constraint for friction cone
                for j=1:obj.nC
                    f(j) = mu(j)*z(1+(1+obj.nD)*(j-1)) - ones(1,obj.nD)*z((2:5)+(1+obj.nD)*(j-1));
                    df(j,nq+nv+obj.nC+1+(1+obj.nD)*(j-1)) = mu(j);
                    df(j,nq+nv+obj.nC+(2:5)+(1+obj.nD)*(j-1)) = -ones(1,obj.nD);
                end
                
                f((1+obj.nC):1+obj.nD:end) = phi;
                df(1+obj.nC:1+obj.nD:end,1:nq) = n;
                for j=1:obj.nD
                    f(1+obj.nC+j:1+obj.nD:end) = gamma+D{j}*v;
                    df(1+obj.nC+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                    df(1+obj.nC+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                    df(1+obj.nC+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
                end
                
                %f_numeric = f;
                %df_numeric = df;
                
                % [f_numeric,df_numeric] = geval(@(y) ERM_nonlincompl_fun_check(y),y,struct('grad_method','numerical'));
                % valuecheck(df,df_numeric,1e-5);
                % valuecheck(f,f_numeric,1e-5);
                
                function [f,df] = ERM_nonlincompl_fun_check(y)
                    nq = obj.plant.getNumPositions;
                    nv = obj.plant.getNumVelocities;
                    x = y(1:nq+nv+obj.nC);
                    z = y(nq+nv+obj.nC+1:end);
                    gamma = x(nq+nv+1:end);
                    q = x(1:nq);
                    v = x(nq+1:nq+nv);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        if isempty(terrain_index)
                            [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                        else
                            [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                        end
                    elseif isempty(manip.uncertainty_source)
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    end
                    
                    %% debugging
                    % a test of using q0 to derive v1, instead of using v1
                    % x0 = zeros(12,1);
                    % q0 = x0(1:nq);%zeros(12,1);
                    % v0 = x0(nq+1:nq+nv);
                    % h = 0.01;
                    % u = zeros(3,1);
                    % [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    %
                    % [M,C,B,dM,dC,dB] = obj.plant.manipulatorDynamics(q,v);
                    % Minv = inv(M);
                    %
                    %v1_est = v0 + Minv*(B*u - C + D{1}'*[z(2);z(5)] + D{2}'*[z(3);z(6)] + n'*[z(1);z(4)])*h;
                    %% end debugging
                    
                    f = zeros(obj.nC*(2+obj.nD),1);
                    df = zeros(obj.nC*(2+obj.nD),nq+nv+obj.nC*(2+obj.nD));
                    
                    %add linear complementarity constraint for friction cone
                    for j=1:obj.nC
                        f(j) = mu(j)*z(1+(1+obj.nD)*(j-1)) - ones(1,obj.nD)*z((2:5)+(1+obj.nD)*(j-1));
                        df(j,nq+nv+obj.nC+1+(1+obj.nD)*(j-1)) = mu(j);
                        df(j,nq+nv+obj.nC+(2:5)+(1+obj.nD)*(j-1)) = -ones(1,obj.nD);
                    end
                    
                    f((1+obj.nC):1+obj.nD:end) = phi;
                    df(1+obj.nC:1+obj.nD:end,1:nq) = n;
                    for j=1:obj.nD
                        f(1+obj.nC+j:1+obj.nD:end) = gamma+D{j}*v;
                        df(1+obj.nC+j:1+obj.nD:end,nq+nv+(1:obj.nC)) = eye(size(D{j},1));  %d/dgamma
                        df(1+obj.nC+j:1+obj.nD:end,nq+(1:nv)) = D{j};%d/dv
                        df(1+obj.nC+j:1+obj.nD:end,1:nq) = matGradMult(dD{j},v);%d/dq
                    end
                end
                
                persistent LCP_non_robust_NCP_residual
                if obj.verbose_print == 1
                    NCP_residual = f.*z;
                    % compute the sum of tangential components
                    NCP_residual_tangential = sum(NCP_residual(2:3));
                    NCP_residual_tangential = NCP_residual_tangential + sum(NCP_residual(5:6));
                    %disp('non-robust NCP residual square');
                    LCP_non_robust_NCP_residual = [LCP_non_robust_NCP_residual NCP_residual_tangential^2];
                    if length(LCP_non_robust_NCP_residual) == obj.N-1
                        %LCP_non_robust_NCP_residual
                        LCP_non_robust_NCP_residual = [];
                    end
                end
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
                
                if strcmp(manip.uncertainty_source, 'friction_coeff')
                    uncertain_mu = [];%reset mu, make it average
                    [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                elseif strcmp(manip.uncertainty_source, 'terrain_height')
                    sample_num = length(manip.uncertain_position_set);
                    for i=1:sample_num
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        terrain_index = i;
                        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                        phi_stack(:,i) = phi;
                        n_stack(:,:,i) = n;
                        dn_stack(:,:,i) = dn;
                        D_stack1(:,:,i) = D{1};
                        D_stack2(:,:,i) = D{2};
                        D_stack3(:,:,i) = D{3};
                        D_stack4(:,:,i) = D{4};
                        dD_stack1(:,:,i) = dD{1};
                        dD_stack2(:,:,i) = dD{2};
                        dD_stack3(:,:,i) = dD{3};
                        dD_stack4(:,:,i) = dD{4};
                    end
                    phi = sum(phi_stack,2)/sample_num;
                    n = sum(n_stack,3)/sample_num;
                    dn = sum(dn_stack,3)/sample_num;
                    D{1} = sum(D_stack1,3)/sample_num;
                    D{2} = sum(D_stack2,3)/sample_num;
                    D{3} = sum(D_stack3,3)/sample_num;
                    D{4} = sum(D_stack4,3)/sample_num;
                    dD{1} = sum(dD_stack1,3)/sample_num;
                    dD{2} = sum(dD_stack2,3)/sample_num;
                    dD{3} = sum(dD_stack3,3)/sample_num;
                    dD{4} = sum(dD_stack4,3)/sample_num;
                elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                    sample_num = length(manip.uncertain_mu_set);
                    for i=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                        manip.uncertain_phi = manip.uncertain_position_set(:,i);
                        uncertain_mu = manip.uncertain_mu_set(i);
                        terrain_index = i;
                        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                        phi_stack(:,i) = phi;
                        n_stack(:,:,i) = n;
                        dn_stack(:,:,i) = dn;
                        D_stack1(:,:,i) = D{1};
                        D_stack2(:,:,i) = D{2};
                        D_stack3(:,:,i) = D{3};
                        D_stack4(:,:,i) = D{4};
                        dD_stack1(:,:,i) = dD{1};
                        dD_stack2(:,:,i) = dD{2};
                        dD_stack3(:,:,i) = dD{3};
                        dD_stack4(:,:,i) = dD{4};
                    end
                    phi = sum(phi_stack,2)/sample_num;
                    n = sum(n_stack,3)/sample_num;
                    dn = sum(dn_stack,3)/sample_num;
                    D{1} = sum(D_stack1,3)/sample_num;
                    D{2} = sum(D_stack2,3)/sample_num;
                    D{3} = sum(D_stack3,3)/sample_num;
                    D{4} = sum(D_stack4,3)/sample_num;
                    dD{1} = sum(dD_stack1,3)/sample_num;
                    dD{2} = sum(dD_stack2,3)/sample_num;
                    dD{3} = sum(dD_stack3,3)/sample_num;
                    dD{4} = sum(dD_stack4,3)/sample_num;
                elseif isempty(manip.uncertainty_source)%non-robust version
                    [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                else
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                    save -ascii friction_coeff_noise.dat w_mu
                end
                
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
                
                %f_numeric = f;
                %df_numeric = df;
                
                % [f_numeric,df_numeric] = geval(@(y) nonlincompl_fun_check(y),y,struct('grad_method','numerical'));
                % valuecheck(df,df_numeric,1e-5);
                % valuecheck(f,f_numeric,1e-5);
                
                function [f,df] = nonlincompl_fun_check(y)
                    nq = obj.plant.getNumPositions;
                    nv = obj.plant.getNumVelocities;
                    x = y(1:nq+nv+obj.nC);
                    z = y(nq+nv+obj.nC+1:end);
                    gamma = x(nq+nv+1:end);
                    
                    q = x(1:nq);
                    v = x(nq+1:nq+nv);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = [];%reset mu, make it average
                        [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        sample_num = length(manip.uncertain_position_set);
                        for i=1:sample_num
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            terrain_index = i;
                            [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                            phi_stack(:,i) = phi;
                            n_stack(:,:,i) = n;
                            dn_stack(:,:,i) = dn;
                            D_stack1(:,:,i) = D{1};
                            D_stack2(:,:,i) = D{2};
                            D_stack3(:,:,i) = D{3};
                            D_stack4(:,:,i) = D{4};
                            dD_stack1(:,:,i) = dD{1};
                            dD_stack2(:,:,i) = dD{2};
                            dD_stack3(:,:,i) = dD{3};
                            dD_stack4(:,:,i) = dD{4};
                        end
                        phi = sum(phi_stack,2)/sample_num;
                        n = sum(n_stack,3)/sample_num;
                        dn = sum(dn_stack,3)/sample_num;
                        D{1} = sum(D_stack1,3)/sample_num;
                        D{2} = sum(D_stack2,3)/sample_num;
                        D{3} = sum(D_stack3,3)/sample_num;
                        D{4} = sum(D_stack4,3)/sample_num;
                        dD{1} = sum(dD_stack1,3)/sample_num;
                        dD{2} = sum(dD_stack2,3)/sample_num;
                        dD{3} = sum(dD_stack3,3)/sample_num;
                        dD{4} = sum(dD_stack4,3)/sample_num;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        sample_num = length(manip.uncertain_mu_set);
                        for i=1:sample_num
                            manip.uncertain_mu = manip.uncertain_mu_set(i);% uncertainy in mu is streamed down into the contactConstraint_manual function
                            manip.uncertain_phi = manip.uncertain_position_set(:,i);
                            uncertain_mu = manip.uncertain_mu_set(i);
                            terrain_index = i;
                            [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q,false,obj.options.active_collision_options);
                            phi_stack(:,i) = phi;
                            n_stack(:,:,i) = n;
                            D_stack1(:,:,i) = D{1};
                            D_stack2(:,:,i) = D{2};
                            D_stack3(:,:,i) = D{3};
                            D_stack4(:,:,i) = D{4};
                            dD_stack1(:,:,i) = dD{1};
                            dD_stack2(:,:,i) = dD{2};
                            dD_stack3(:,:,i) = dD{3};
                            dD_stack4(:,:,i) = dD{4};
                        end
                        phi = sum(phi_stack,2)/sample_num;
                        n = sum(n_stack,3)/sample_num;
                        dn = sum(dn_stack,3)/sample_num;
                        D{1} = sum(D_stack1,3)/sample_num;
                        D{2} = sum(D_stack2,3)/sample_num;
                        D{3} = sum(D_stack3,3)/sample_num;
                        D{4} = sum(D_stack4,3)/sample_num;
                        dD{1} = sum(dD_stack1,3)/sample_num;
                        dD{2} = sum(dD_stack2,3)/sample_num;
                        dD{3} = sum(dD_stack3,3)/sample_num;
                        dD{4} = sum(dD_stack4,3)/sample_num;
                    elseif isempty(manip.uncertainty_source)%non-robust version
                        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q,false,obj.options.active_collision_options);
                    else
                        w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);%friction coefficient noise
                        save -ascii friction_coeff_noise.dat w_mu
                    end
                    
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
        end
        
        function [c,dc] = robustVariancecost_ML(obj, x_full, Fext_full)
            global timestep_updated
            global x_previous
            global df_previous
            global terrain_index
            global uncertain_mu
            manip = obj.plant.getManipulator();
            
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.nFext;%used to be obj.plant.getNumInputs;
            obj.nu = nu;
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px_init = obj.cached_Px(:,:,1);
            
            if strcmp(manip.uncertainty_source,'friction_coeff')
                w_mu = manip.uncertain_mu_set;
                w_noise = [w_mu];
                Pw = diag([0.01]);
            elseif strcmp(manip.uncertainty_source,'terrain_height')
                w_phi = manip.uncertain_position_set;
                w_noise = [w_phi];
                Pw = diag([0.016]);
            elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                w_mu = manip.uncertain_mu_set;
                w_phi = manip.uncertain_position_set;
                w_noise = [w_mu;w_phi];
                Pw = diag([0.01, 0.016]);
            elseif isempty(manip.uncertainty_source)
                Pw = [];
                w_noise = [];
            elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                save -ascii friction_coeff_noise.dat w_mu
                %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                %w_phi = [x;y];%height noise
                %save -ascii initial_position_noise.dat w_phi
                w_noise = [w_mu];
            end
            
            % disturbance variance
            % currently only consider terrain height and/or friction coefficient uncertainty
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            nw = size(Pw,1);
            sampling_method = 2;%option 1: unscented transform, option 2: random sampling with a smaller number
            if sampling_method == 1
                n_sampling_point = 2*(obj.nx+nw);
            elseif sampling_method == 2
                n_sampling_point = 15;
                w_state = load('state_noise.dat');
                %w_state = 0.001*randn(12,50);
                %save -ascii state_noise.dat w_state
            end
            
            w_avg = 1/n_sampling_point;
            K = obj.options.K;%controller gain
            
            %initialize cost c and dc
            kappa = obj.options.kappa;%cost scaling parameter
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix c(k=1) = Px(1);
            %c = 0;
            c = kappa*trace(obj.cached_Px(:,:,1));
            %c = 1/2*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(12)));
            dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(timestep_updated,Sig,u_fdb_k)
                [xdn,df] = obj.plant.update(timestep_updated,Sig,u_fdb_k);
            end
            
            plant_update = @objPlantUpdate;
            %obj.plant.time_step = time_step;
            %option 1: one sigma has one environment sample; 
            %option 2: one sigma runs through all environment samples (full version, computationally heavy)
            noise_sample_type = 1;
            df_previous_full = [];
            xdn_previous_full = [];
             
            for k = 1:obj.N-1
                if sampling_method == 1
                    %covariance matrix set-up
                    %always assume constant covariance matrix, there is another
                    %cost version, i.e., robustVariancecost() uses the
                    %covariance computed from previous time step. But more
                    %chain rule components have to be added.
                    Px_init = obj.cached_Px(:,:,1);
                    [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
                    
                    if d
                        diverge = k;
                        return;
                    end
                    S = scale*S;
                    Sig_init(:,:,k) = [S -S];
                    if k == 1%customize for ML formulation
                        Sig(:,:,k) = Sig_init(:,:,k);
                    end
                    
                    %state dimensionality of sigma point depends on the number of uncertainties.
                    %Our case w_noise has 0, 1, or 2 dimensions
                    if isempty(w_noise)
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
                        end
                    else
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
                        end
                    end
                    
                    % double check that the mean of sigma points at the first step is zero
                    if k == 1
                        x_mean(:,k) = zeros(obj.nx,1);
                        for j = 1:n_sampling_point
                            x_mean(:,k) = x_mean(:,k) + w_avg*Sig_init(1:obj.nx,j,k);
                        end
                        %c = c + norm(x(:,k)-x_mean(:,k))^2;% this mean is zero for the first time step
                    end
                elseif sampling_method == 2
                    for j = 1:n_sampling_point
                        Sig_init(:,j,k) = x(:,k) + w_state(:,j);
                    end
                end
                 
                %Propagate sigma points through nonlinear dynamics
                for j = 1:n_sampling_point                    
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:obj.nx/2,j,k),Sig_init(obj.nx/2+1:obj.nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = w_mu(j);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        terrain_index = j;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        uncertain_mu = w_mu(j);
                        terrain_index = j;
                    end
                    
                    % add feedback control
                    t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig_init(1:obj.nx,j,k) - x(:,k));
                    
                    if k > 1
                        x_previous = xdn_previous_full(:,j);
                        df_previous = df_previous_full(:,:,j);
                    end
                    
                    %% two different noise sampling methods
                    %option 1: one sigma has one environment sample; 
                    %option 2: one sigma runs through all environment samples (full version, computationally heavy)
                    if noise_sample_type == 1
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                    elseif noise_sample_type == 2
                        xdn_analytical(:,j) = zeros(nx,1);
                        df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                        for kk=1:length(w_noise)
                            [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                            xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                            df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                        end
                        xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                        df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                    end
                    
                    %% numerical diff for debugging
                    % dt = diag(max(sqrt(eps(timestep_updated)),1e-7));
                    % dx = diag(max(sqrt(eps(Sig_init(1:obj.nx,j,k))), 1e-7));
                    % du = diag(max(sqrt(eps(u_fdb_k)),1e-7));
                    %
                    % [xdnp,~] = feval(plant_update,timestep_updated+dt,Sig_init(1:obj.nx,j,k),u_fdb_k);
                    % [xdnm,~] = feval(plant_update,timestep_updated-dt,Sig_init(1:obj.nx,j,k),u_fdb_k);
                    % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                    %
                    % N_finite_diff_x = length(Sig_init(1:obj.nx,j,k));
                    % for m = 1:N_finite_diff_x
                    %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                    %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                    %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                    % end
                    %
                    % N_finite_diff_u = length(u_fdb_k);
                    % for m = 1:N_finite_diff_u
                    %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k),u_fdb_k+du(:,m));
                    %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k),u_fdb_k-du(:,m));
                    %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                    % end
                    
                    % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                    %     keyboard
                    % end
                    %for jjj=1:14
                    %    jjj
                    %    [df_analytical(:,jjj,j),df_numeric(:,jjj)]'
                    %end
                    
                    %% update state and dynamics gradient 
                    xdn = xdn_analytical;
                    df(:,:,j) = df_analytical(:,:,j);
                    
                    Sig(1:nx,j,k+1) = xdn(1:nx,j);
                    dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                    dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                    dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                end
                xdn_previous_full = xdn;
                df_previous_full = df;
                
                %% calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(obj.nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                end
                
                %% shift sigma point
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                %         eta_final(k+1) = eta_optimal;
                %     end
                % end
                %
                % % scale Sigma point deviation by eta
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                % end
                %
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(obj.nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                % end
                %
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                
                Px(:,:,k+1) = zeros(obj.nx);
                for jj = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                end
                 
                %% accumulate returned cost (several different options)
                %c = c + kappa*trace(Px(:,:,k+1));
                %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12)));%ML mean deviation version
                c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                
                %% derivative of variance matrix
                % gradient of MLcost w.r.t state x and control u
                dMLcostdx(:,:,k+1) = zeros(obj.N,obj.nx);
                dMLcostdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                % sum of sigma mean gradient w.r.t state x
                dSig_m_kplus1_dx_sum = zeros(obj.nx);
                % sum of sigma mean gradient w.r.t control u
                dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                dSig_i_kplus1_dx_set = zeros(obj.nx, obj.nx, n_sampling_point);
                dSig_i_kplus1_du_set = zeros(obj.nx, nu, n_sampling_point);
                
                for i=1:n_sampling_point
                    if i == 1
                        % this for-loop is for sigma mean \mu{x_tilde}_{k+1}, 
                        % and only needs to go through once since the mean remains the same for different sigma points
                        for m = 1:n_sampling_point
                            % sum of sigma gradient w.r.t control x and u
                            dSig_m_kplus1_dx = zeros(obj.nx);
                            dSig_m_kplus1_du = zeros(obj.nx,1);
                            
                            dSig_m_kplus1_dx = dfdx(:,:,m,k+1) + dfdSig(:,:,m,k+1);
                            dSig_m_kplus1_du = dfdu(:,:,m,k+1);% [double check that du is not affected]
                            
                            dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                            dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                        end
                    end
                    
                    % run n_sampling_point times in total to obtain gradient w.r.t sigma points
                    dSig_i_kplus1_dx = zeros(obj.nx);
                    dSig_i_kplus1_du = zeros(obj.nx,nu);
                    
                    dSig_i_kplus1_dx = dfdx(:,:,i,k+1) + dfdSig(:,:,i,k+1);
                    dSig_i_kplus1_du = dfdu(:,:,i,k+1);
                    
                    dSig_i_kplus1_dx_set(:,:,i) = dSig_i_kplus1_dx;
                    dSig_i_kplus1_du_set(:,:,i) = dSig_i_kplus1_du;
                end
                
                % one gradient component of first quadratic ML term w.r.t x and u
                dMLcostdx(k,:,k+1) = dMLcostdx(k,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                dMLcostdu(k,:,k+1) = dMLcostdu(k,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                
                % one gradient component of first quadratic ML term w.r.t x and u
                % and gradient component of second quadratic ML term w.r.t x and u
                % this part has to go through element by element.
                dCovdx = zeros(nx,nx,nx);
                dCovdu = zeros(nx,nx,nu);
                
                for pp=1:nx
                    for i=1:n_sampling_point
                        dCovdmeandev = zeros(nx,obj.nx,nx);%d_covariance /d_(sigma mean deviation from nominal state)
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
                            dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_set(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
                            if pp <= nu
                                dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_set(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
                            end
                        end
                    end
                    % gradient w.r.t x of the first quadratic term in ML cost
                    dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdx(:,:,pp));
                    % (1st option) gradient w.r.t x of the second quadratic term in ML cost
                    %dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdx(:,:,pp));
                    % (2nd option) gradient w.r.t x of the second quadratic term in ML cost
                    dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) + trace(dCovdx(:,:,pp));
                    if pp <= nu
                        % gradient w.r.t u of the first quadratic term in ML cost
                        dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdu(:,:,pp));
                        % (1st option) gradient w.r.t u of the second quadratic term in ML cost
                        %dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdu(:,:,pp));
                        % (2nd option) gradient w.r.t u of the second quadratic term in ML cost
                        dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) + trace(dCovdu(:,:,pp));
                    end
                end
            end
             
            %% reaccumulate all gradient components at all N time steps
            dc = [];
            for jj=1:obj.N % for state gradient
                if (jj == 1)
                    dMLcost_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                else
                    %dMLcost_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                    dMLcost_sum_dx_k = (pinv(Px(:,:,jj)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                end
                
                if jj < obj.N
                    dMLcost_sum_dx_k = dMLcost_sum_dx_k + dMLcostdx(jj,:,jj+1);
                end
                dc = [dc, dMLcost_sum_dx_k];
            end
            
            for jj=1:obj.N % for control gradient
                dMLcost_sum_du_k = zeros(1, nu);
                if jj < obj.N
                    dMLcost_sum_du_k = dMLcost_sum_du_k + dMLcostdu(jj,:,jj+1);
                end
                dc = [dc, dMLcost_sum_du_k];
            end
            
            %% scale robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
            fprintf('ML robust cost: %4.4f\n',c);
            tElapsed = toc(tStart);
            
            %% check numerical gradient
            X0 = [x_full; Fext_full];
            c_numeric = c;
            dc_numeric = dc;
            % [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_ML_check(obj,X0),X0,struct('grad_method','numerical'));
            % [c_numeric,dc_numeric] = robustVariancecost_check(obj, X0);
            % valuecheck(dc,dc_numeric,1e-5);
            % valuecheck(c,c_numeric,1e-5);
            
            %% another way to check numerical gradient
            % X0 = X0 + randn(size(X0))*0.1;
            % fun = @(X0) robustVariancecost_check(obj, X0);
            % DerivCheck(fun, X0)
            
            function DerivCheck(funptr, X0, ~, varargin)
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %  Checks the analytic gradient of a function 'funptr' at a point X0, and
                %  compares to numerical gradient.  Useful for checking gradients computed
                %  for fminunc and fmincon.
                %  Call with same arguments as you would call for optimization (fminunc).
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
                dd =  dot(rr, JJ)-f2+f1;
                fprintf('difference between numerical and analytical: %4.15f\n',dd);
            end
            
            %% for numerical gradient check
            function [c,dc] = robustVariancecost_ML_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                manip = obj.plant.getManipulator();
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nx = nq+nv;
                nu = obj.nFext;%used to be obj.plant.getNumInputs;
                obj.nu = nu;
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px_init = obj.cached_Px(:,:,1);
                
                if strcmp(manip.uncertainty_source,'friction_coeff')
                    w_mu = manip.uncertain_mu_set;
                    w_noise = [w_mu];
                    Pw = diag([0.01]);
                elseif strcmp(manip.uncertainty_source,'terrain_height')
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_phi];
                    Pw = diag([0.016]);
                elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                    w_mu = manip.uncertain_mu_set;
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_mu;w_phi];
                    Pw = diag([0.01, 0.016]);
                elseif isempty(manip.uncertainty_source)
                    Pw = [];
                    w_noise = [];
                elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                    %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                    w_noise = [w_mu];
                end
                
                % disturbance variance
                % currently only consider terrain height and/or friction coefficient uncertainty
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
                nw = size(Pw,1);
                sampling_method = 2;%option 1: unscented transform, option 2: random sampling with a smaller number
                if sampling_method == 1
                    n_sampling_point = 2*(obj.nx+nw);
                elseif sampling_method == 2
                    n_sampling_point = 6;
                end
                
                w_avg = 1/n_sampling_point;
                K = obj.options.K;%controller gain
                
                %initialize cost c and dc
                kappa = obj.options.kappa;%cost scaling parameter
                x_mean = zeros(obj.nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix c(k=1) = Px(1);
                %c = 0;
                c = kappa*trace(obj.cached_Px(:,:,1));
                %c = 1/2*log(det(Px(:,:,1)+obj.options.Px_regularizer_coeff*eye(12)));
                dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(timestep_updated,Sig,u_fdb_k)
                    [xdn,df] = obj.plant.update(timestep_updated,Sig,u_fdb_k);
                end
                
                plant_update = @objPlantUpdate;
                %obj.plant.time_step = time_step;
                %option 1: one sigma has one environment sample;
                %option 2: one sigma runs through all environment samples (full version, computationally heavy)
                noise_sample_type = 1;
                df_previous_full = [];
                xdn_previous_full = [];
                
                for k = 1:obj.N-1
                    if sampling_method == 1
                        %covariance matrix set-up
                        %always assume constant covariance matrix, there is another
                        %cost version, i.e., robustVariancecost() uses the
                        %covariance computed from previous time step. But more
                        %chain rule components have to be added.
                        Px_init = obj.cached_Px(:,:,1);
                        [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
                        
                        if d
                            diverge = k;
                            return;
                        end
                        S = scale*S;
                        Sig_init(:,:,k) = [S -S];
                        if k == 1%customize for ML formulation
                            Sig(:,:,k) = Sig_init(:,:,k);
                        end
                        
                        %state dimensionality of sigma point depends on the number of uncertainties.
                        %Our case w_noise has 0, 1, or 2 dimensions
                        if isempty(w_noise)
                            for j = 1:n_sampling_point
                                Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
                            end
                        else
                            for j = 1:n_sampling_point
                                Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
                            end
                        end
                        
                        % double check that the mean of sigma points at the first step is zero
                        if k == 1
                            x_mean(:,k) = zeros(obj.nx,1);
                            for j = 1:n_sampling_point
                                x_mean(:,k) = x_mean(:,k) + w_avg*Sig_init(1:obj.nx,j,k);
                            end
                            %c = c + norm(x(:,k)-x_mean(:,k))^2;% this mean is zero for the first time step
                        end
                    elseif sampling_method == 2
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = x(:,k);% + 0.001*randn(nx,1);
                        end
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:n_sampling_point
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:obj.nx/2,j,k),Sig_init(obj.nx/2+1:obj.nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        
                        if strcmp(manip.uncertainty_source, 'friction_coeff')
                            uncertain_mu = w_mu(j);
                        elseif strcmp(manip.uncertainty_source, 'terrain_height')
                            terrain_index = j;
                        elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                            uncertain_mu = w_mu(j);
                            terrain_index = j;
                        end
                        
                        % add feedback control
                        t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                        u_fdb_k = u(:,k) - K*(Sig_init(1:obj.nx,j,k) - x(:,k));
                        
                        if k > 1
                            x_previous = xdn_previous_full(:,j);
                            df_previous = df_previous_full(:,:,j);
                        end
                        
                        %% two different noise sampling methods
                        %option 1: one sigma has one environment sample;
                        %option 2: one sigma runs through all environment samples (full version, computationally heavy)
                        if noise_sample_type == 1
                            [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                        elseif noise_sample_type == 2
                            xdn_analytical(:,j) = zeros(nx,1);
                            df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                            for kk=1:length(w_noise)
                                [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                                xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                                df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                            end
                            xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                            df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                        end
                        
                        %% numerical diff for debugging
                        % dt = diag(max(sqrt(eps(timestep_updated)),1e-7));
                        % dx = diag(max(sqrt(eps(Sig_init(1:obj.nx,j,k))), 1e-7));
                        % du = diag(max(sqrt(eps(u_fdb_k)),1e-7));
                        %
                        % [xdnp,~] = feval(plant_update,timestep_updated+dt,Sig_init(1:obj.nx,j,k),u_fdb_k);
                        % [xdnm,~] = feval(plant_update,timestep_updated-dt,Sig_init(1:obj.nx,j,k),u_fdb_k);
                        % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                        %
                        % N_finite_diff_x = length(Sig_init(1:obj.nx,j,k));
                        % for m = 1:N_finite_diff_x
                        %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                        %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                        %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                        % end
                        %
                        % N_finite_diff_u = length(u_fdb_k);
                        % for m = 1:N_finite_diff_u
                        %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k),u_fdb_k+du(:,m));
                        %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:obj.nx,j,k),u_fdb_k-du(:,m));
                        %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                        % end
                        
                        % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                        %     keyboard
                        % end
                        %for jjj=1:14
                        %    jjj
                        %    [df_analytical(:,jjj,j),df_numeric(:,jjj)]'
                        %end
                        
                        %% update state and dynamics gradient
                        xdn = xdn_analytical;
                        df(:,:,j) = df_analytical(:,:,j);
                        
                        Sig(1:nx,j,k+1) = xdn(1:nx,j);
                        dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                        dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                    end
                    xdn_previous_full = xdn;
                    df_previous_full = df;
                    
                    %% calculate mean and variance w.r.t. [x_k] from sigma points
                    x_mean(:,k+1) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    end
                    
                    %% shift sigma point
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                    %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                    %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                    %         eta_final(k+1) = eta_optimal;
                    %     end
                    % end
                    %
                    % % scale Sigma point deviation by eta
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                    % end
                    %
                    % % recalculate mean and variance w.r.t. [x_k] from sigma points
                    % x_mean(:,k+1) = zeros(obj.nx,1);
                    % for j = 1:n_sampling_point
                    %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    % end
                    %
                    % % check that the mean deviation term is cancelled out
                    % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                    %     disp('shifting scheme is not correct')
                    %     keyboard
                    % end
                    
                    Px(:,:,k+1) = zeros(obj.nx);
                    for jj = 1:n_sampling_point
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    end
                    
                    %% accumulate returned cost (several different options)
                    %c = c + kappa*trace(Px(:,:,k+1));
                    %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12)));%ML mean deviation version
                    c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)) + trace(Px(:,:,k+1));%ML mean deviation version
                    
                    %% derivative of variance matrix
                    % gradient of MLcost w.r.t state x and control u
                    dMLcostdx(:,:,k+1) = zeros(obj.N,obj.nx);
                    dMLcostdu(:,:,k+1) = zeros(obj.N-1,nu);
                    
                    % sum of sigma mean gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(obj.nx);
                    % sum of sigma mean gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                    dSig_i_kplus1_dx_set = zeros(obj.nx, obj.nx, n_sampling_point);
                    dSig_i_kplus1_du_set = zeros(obj.nx, nu, n_sampling_point);
                    
                    for i=1:n_sampling_point
                        if i == 1
                            % this for-loop is for sigma mean \mu{x_tilde}_{k+1},
                            % and only needs to go through once since the mean remains the same for different sigma points
                            for m = 1:n_sampling_point
                                % sum of sigma gradient w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(obj.nx);
                                dSig_m_kplus1_du = zeros(obj.nx,1);
                                
                                dSig_m_kplus1_dx = dfdx(:,:,m,k+1) + dfdSig(:,:,m,k+1);
                                dSig_m_kplus1_du = dfdu(:,:,m,k+1);% [double check that du is not affected]
                                
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run n_sampling_point times in total to obtain gradient w.r.t sigma points
                        dSig_i_kplus1_dx = zeros(obj.nx);
                        dSig_i_kplus1_du = zeros(obj.nx,nu);
                        
                        dSig_i_kplus1_dx = dfdx(:,:,i,k+1) + dfdSig(:,:,i,k+1);
                        dSig_i_kplus1_du = dfdu(:,:,i,k+1);
                        
                        dSig_i_kplus1_dx_set(:,:,i) = dSig_i_kplus1_dx;
                        dSig_i_kplus1_du_set(:,:,i) = dSig_i_kplus1_du;
                    end
                    
                    % one gradient component of first quadratic ML term w.r.t x and u
                    dMLcostdx(k,:,k+1) = dMLcostdx(k,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    dMLcostdu(k,:,k+1) = dMLcostdu(k,:,k+1) + (pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    
                    % one gradient component of first quadratic ML term w.r.t x and u
                    % and gradient component of second quadratic ML term w.r.t x and u
                    % this part has to go through element by element.
                    dCovdx = zeros(nx,nx,nx);
                    dCovdu = zeros(nx,nx,nu);
                    
                    for pp=1:nx
                        for i=1:n_sampling_point
                            dCovdmeandev = zeros(nx,obj.nx,nx);%d_covariance /d_(sigma mean deviation from nominal state)
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
                                dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_set(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
                                if pp <= nu
                                    dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_set(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
                                end
                            end
                        end
                        % gradient w.r.t x of the first quadratic term in ML cost
                        dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdx(:,:,pp));
                        % (1st option) gradient w.r.t x of the second quadratic term in ML cost
                        %dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdx(:,:,pp));
                        % (2nd option) gradient w.r.t x of the second quadratic term in ML cost
                        dMLcostdx(k,pp,k+1) = dMLcostdx(k,pp,k+1) + trace(dCovdx(:,:,pp));
                        if pp <= nu
                            % gradient w.r.t u of the first quadratic term in ML cost
                            dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdu(:,:,pp));
                            % (1st option) gradient w.r.t u of the second quadratic term in ML cost
                            %dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1)+obj.options.Px_regularizer_coeff*eye(12))*dCovdu(:,:,pp));
                            % (2nd option) gradient w.r.t u of the second quadratic term in ML cost
                            dMLcostdu(k,pp,k+1) = dMLcostdu(k,pp,k+1) + trace(dCovdu(:,:,pp));
                        end
                    end
                end
                
                %% reaccumulate all gradient components at all N time steps
                dc = [];
                for jj=1:obj.N % for state gradient
                    if (jj == 1)
                        dMLcost_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                    else
                        %dMLcost_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                        dMLcost_sum_dx_k = (pinv(Px(:,:,jj)+obj.options.Px_regularizer_coeff*eye(12))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                    end
                    
                    if jj < obj.N
                        dMLcost_sum_dx_k = dMLcost_sum_dx_k + dMLcostdx(jj,:,jj+1);
                    end
                    dc = [dc, dMLcost_sum_dx_k];
                end
                
                for jj=1:obj.N % for control gradient
                    dMLcost_sum_du_k = zeros(1, nu);
                    if jj < obj.N
                        dMLcost_sum_du_k = dMLcost_sum_du_k + dMLcostdu(jj,:,jj+1);
                    end
                    dc = [dc, dMLcost_sum_du_k];
                end
                
                %% scale robust cost
                c = obj.options.contact_robust_cost_coeff*c;
                dc = obj.options.contact_robust_cost_coeff*dc;
            end
        end
        
        
        function [c,dc] = robustVariancecost_ML_one_step_prop(obj, x_full, Fext_full)
            global timestep_updated
            global time_step
            global x_previous
            global df_previous
            global terrain_index
            global uncertain_mu
            manip = obj.plant.getManipulator();
                        
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.nFext;%obj.plant.getNumInputs;
            obj.nu = nu;
            u
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px_init = obj.cached_Px(:,:,1);
            
            if strcmp(manip.uncertainty_source,'friction_coeff')
                w_mu = manip.uncertain_mu_set;
                w_noise = [w_mu];
                Pw = diag([0.01]);
            elseif strcmp(manip.uncertainty_source,'terrain_height')
                w_phi = manip.uncertain_position_set;
                w_noise = [w_phi];
                Pw = diag([0.016]);
            elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                w_mu = manip.uncertain_mu_set;
                w_phi = manip.uncertain_position_set;
                w_noise = [w_mu;w_phi];
                Pw = diag([0.01, 0.016]);
            elseif isempty(manip.uncertainty_source)
                Pw = [];
                w_noise = [];
            elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                save -ascii friction_coeff_noise.dat w_mu
                %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                %w_phi = [x;y];%height noise
                %save -ascii initial_position_noise.dat w_phi
                w_noise = [w_mu];
            end
            
            % disturbance variance
            % currently only consider terrain height and/or friction coefficient
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            nw = size(Pw,1);
            n_sampling_point = 1;%2*(obj.nx+nw);
            w_avg = 1/n_sampling_point;
            K = obj.options.K;
            
            %initialize c and dc
            kappa = obj.options.kappa;
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            c = 0;
            %c = kappa*trace(Px(:,:,1));
            dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
            end
            
            plant_update = @objPlantUpdate;
            %obj.plant.time_step = time_step;
            df_previous_full = [];
            xdn_previous_full = [];
            
            for k = 1:obj.N-1%[Ye: double check the index]
                %Propagate sigma points through nonlinear dynamics
                Px_init = obj.cached_Px(:,:,1);%Always assume constant covariance matrix
%                 [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
% 
%                 if d
%                     diverge = k;
%                     return;
%                 end
%                 S = scale*S;
%                 Sig_init(:,:,k) = [S -S];
%                 if k == 1%customize for ML formulation
%                     Sig(:,:,k) = Sig_init(:,:,k);
%                 end
%                 
%                 if isempty(w_noise)
%                     for j = 1:n_sampling_point
%                         Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
%                     end
%                 else
%                     for j = 1:n_sampling_point
%                         Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
%                     end
%                 end
                
                Sig_init(:,1,k) = x(:,k);
                
                if k == 1
                    x_mean(:,k) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k) = x_mean(:,k) + w_avg*Sig_init(1:obj.nx,j,k);
                    end
                    %c = c + norm(x(:,k)-x_mean(:,k))^2;
                end
                
                %Propagate sigma points through nonlinear dynamics
                for j = 1:n_sampling_point
                    noise_index = j;
                    
                    % a hacky way to implement the control input
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:obj.nx/2,j,k),Sig_init(obj.nx/2+1:obj.nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = w_mu(j);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        terrain_index = j;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        uncertain_mu = w_mu(j);
                        terrain_index = j;
                    end
                    
                    % add feedback control
                    t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig_init(1:obj.nx,j,k) - x(:,k));
                    
                    if k > 1
                        x_previous = xdn_previous_full(:,j);
                        df_previous = df_previous_full(:,:,j);
                    end
                    
                    % estimate whether current state is close to contact
                    [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig_init(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                    phi_bottom = phi_current(2:2:end);
                    active_threshold = 0.5;
                    if 0 %any(phi_bottom<active_threshold)
                        for kk = 1:size_terrain_sample
                            if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                obj.plant.friction_coeff = w_mu(kk);
                            elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                obj.plant.terrain_index = kk;
                            end
                            [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                        end
                        xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                        Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                        Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                    else
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig_init(1:nx,j,k),u_fdb_k);
                        
                        % %numerical diff
                        % dt = diag(sqrt(eps(t)));
                        % dx = diag(max(sqrt(eps(Sig(1:obj.nx,j,k))), 1e-7*ones(obj.nx,1)));
                        % du = diag(max(sqrt(eps(u_fdb_k)), 1e-7*ones(nu,1)));
                        %
                        % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                        %
                        % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                        % for m = 1:N_finite_diff_x
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                        %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                        % end
                        %
                        % N_finite_diff_u = length(u_fdb_k);
                        % for m = 1:N_finite_diff_u
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                        %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                        % end
                        %
                        % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                        %     keyboard
                        % end
                        
                        xdn = xdn_analytical;
                        df(:,:,j) = df_analytical(:,:,j);
                        
                        Sig(1:nx,j,k+1) = xdn(1:nx,j);
                        dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                        dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                    end
                end
                xdn_previous_full = xdn;
                df_previous_full = df;
                
                %calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(obj.nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                end
                
                % % shift sigma point
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                %         eta_final(k+1) = eta_optimal;
                %     end
                % end
                %
                % % scale Sigma point deviation by eta
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                % end
                %
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(obj.nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                % end
                %
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                
                Px(:,:,k+1) = zeros(obj.nx);
                for jj = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                end
                %c = c + kappa*trace(Px(:,:,k+1));
                
                % for jj = 1:n_sampling_point
                % V_comp(:,:,k+1) = (Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                % c = c + kappa*trace(w*V_comp(:,:,k+1));
                % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                
                % accumulate returned cost
                %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*(x(:,k+1)-x_mean(:,k+1));%ML mean deviation version
                
                % derivative of variance matrix
                % gradient of Tr(V) w.r.t state vector x
                dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                tic
                % gradient w.r.t state x
                dSig_m_kplus1_dx_sum = zeros(obj.nx);
                % gradient w.r.t control u
                dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                dSig_i_kplus1_dx_set = zeros(obj.nx, obj.nx, n_sampling_point);
                dSig_i_kplus1_du_set = zeros(obj.nx, nu, n_sampling_point);
                
                for i=1:n_sampling_point
                    if i == 1
                        for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                            % gradient of Tr(V_{k+1}) w.r.t control x and u
                            dSig_m_kplus1_dx = zeros(obj.nx);
                            dSig_m_kplus1_du = zeros(obj.nx,1);
                            
                            dSig_m_kplus1_dx = dfdx(:,:,m,k+1) + dfdSig(:,:,m,k+1);
                            dSig_m_kplus1_du = dfdu(:,:,m,k+1);% [double check that du is not affected]
                            
                            dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                            dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                        end
                    end
                    
                    % run 2*(obj.nx+nw) times in total to obtain
                    % gradient w.r.t sigma points
                    dSig_i_kplus1_dx = zeros(obj.nx);
                    dSig_i_kplus1_du = zeros(obj.nx,nu);
                    
                    dSig_i_kplus1_dx = dfdx(:,:,i,k+1) + dfdSig(:,:,i,k+1);
                    dSig_i_kplus1_du = dfdu(:,:,i,k+1);
                    
                    dSig_i_kplus1_dx_set(:,:,i) = dSig_i_kplus1_dx;
                    dSig_i_kplus1_du_set(:,:,i) = dSig_i_kplus1_du;
                end
                
                % gradient of mean residual w.r.t state x and control u, assume norm 2
                dmeanRdx(k,:,k+1) = dmeanRdx(k,:,k+1) + ((x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                dmeanRdu(k,:,k+1) = dmeanRdu(k,:,k+1) + ((x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                
                dCovdx = zeros(nx,nx,nx);
                dCovdu = zeros(nx,nx,nu);
                
%                 for pp=1:nx
%                     for i=1:n_sampling_point
%                         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
%                         for nn=1:nx
%                             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
%                             if nn == 1
%                                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
%                                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
%                             elseif nn == nx
%                                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
%                                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
%                             else
%                                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
%                                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
%                                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
%                             end
%                             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_set(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
%                             if pp <= nu
%                                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_set(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
%                             end
%                         end
%                     end
%                     dmeanRdx(k,pp,k+1) = dmeanRdx(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
%                     dmeanRdx(k,pp,k+1) = dmeanRdx(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
%                     if pp <= nu
%                         dmeanRdu(k,pp,k+1) = dmeanRdu(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
%                         dmeanRdu(k,pp,k+1) = dmeanRdu(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
%                     end
%                 end
            end
            
            dc = [];
            % cost gradient w.r.t x at first time step is zero
            for jj=1:obj.N % index for x_k
                if (jj == 1)
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                else
                    %dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                    dmeanR_sum_dx_k = ((x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                end
                
                if jj < obj.N
                    dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jj+1);
                end
                %dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                dc = [dc, dmeanR_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            for jj=1:obj.N % index for u_k
                dmeanR_sum_du_k = zeros(1, nu);
                if jj < obj.N
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jj+1);
                end
                %dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                dc = [dc, dmeanR_sum_du_k];
            end
            
            % scale this robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
            fprintf('ML robust cost: %4.4f\n',c);
            tElapsed = toc(tStart);
            
            % % check gradient
            % disp('check gradient')
            % c_numeric = c;
            % dc_numeric = dc;
            %
            X0 = [x_full; Fext_full];
            % %X0 = X0 + randn(size(X0))*0.1;
            %
            % fun = @(X0) robustVariancecost_check(obj, X0);
            % DerivCheck(fun, X0)
            % disp('finish')
            c_numeric = c;
            dc_numeric = dc;
            
            % [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_ML_one_step_prop_check(obj,X0),X0,struct('grad_method','numerical'));
            %
            % [c_numeric,dc_numeric] = robustVariancecost_check(obj, X0);
            %
            % valuecheck(dc,dc_numeric,1e-5);
            % valuecheck(c,c_numeric,1e-5);
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %
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
                dd =  dot(rr, JJ)-f2+f1;
                fprintf('difference between numerical and analytical: %4.15f\n',dd);
            end
            
            function [c,dc] = robustVariancecost_ML_one_step_prop_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                
                manip = obj.plant.getManipulator();
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nx = nq+nv;
                nu = obj.nFext;%obj.plant.getNumInputs;
                obj.nu = nu;
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px_init = obj.cached_Px(:,:,1);
                
                if strcmp(manip.uncertainty_source,'friction_coeff')
                    w_mu = manip.uncertain_mu_set;
                    w_noise = [w_mu];
                    Pw = diag([0.01]);
                elseif strcmp(manip.uncertainty_source,'terrain_height')
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_phi];
                    Pw = diag([0.016]);
                elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                    w_mu = manip.uncertain_mu_set;
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_mu;w_phi];
                    Pw = diag([0.01, 0.016]);
                elseif isempty(manip.uncertainty_source)
                    Pw = [];
                    w_noise = [];
                elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                    %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                    w_noise = [w_mu];
                end
                
                % disturbance variance
                % currently only consider terrain height and/or friction coefficient
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
                nw = size(Pw,1);
                n_sampling_point = 1;%2*(obj.nx+nw);
                w_avg = 1/n_sampling_point;
                K = obj.options.K;
                
                %initialize c and dc
                kappa = obj.options.kappa;
                x_mean = zeros(obj.nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
                c = 0;
                %c = kappa*trace(Px(:,:,1));
                dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                    [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
                end
                
                plant_update = @objPlantUpdate;
                %obj.plant.time_step = time_step;
                df_previous_full = [];
                xdn_previous_full = [];
                
                for k = 1:obj.N-1%[Ye: double check the index]
                    %Propagate sigma points through nonlinear dynamics
                    Px_init = obj.cached_Px(:,:,1);%Always assume constant covariance matrix
                    %                 [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
                    %
                    %                 if d
                    %                     diverge = k;
                    %                     return;
                    %                 end
                    %                 S = scale*S;
                    %                 Sig_init(:,:,k) = [S -S];
                    %                 if k == 1%customize for ML formulation
                    %                     Sig(:,:,k) = Sig_init(:,:,k);
                    %                 end
                    %
                    %                 if isempty(w_noise)
                    %                     for j = 1:n_sampling_point
                    %                         Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
                    %                     end
                    %                 else
                    %                     for j = 1:n_sampling_point
                    %                         Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
                    %                     end
                    %                 end
                    
                    Sig_init(:,1,k) = x(:,k);
                    
                    if k == 1
                        x_mean(:,k) = zeros(obj.nx,1);
                        for j = 1:n_sampling_point
                            x_mean(:,k) = x_mean(:,k) + w_avg*Sig_init(1:obj.nx,j,k);
                        end
                        %c = c + norm(x(:,k)-x_mean(:,k))^2;
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:n_sampling_point
                        noise_index = j;
                        
                        % a hacky way to implement the control input
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:obj.nx/2,j,k),Sig_init(obj.nx/2+1:obj.nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        
                        if strcmp(manip.uncertainty_source, 'friction_coeff')
                            uncertain_mu = w_mu(j);
                        elseif strcmp(manip.uncertainty_source, 'terrain_height')
                            terrain_index = j;
                        elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                            uncertain_mu = w_mu(j);
                            terrain_index = j;
                        end
                        
                        % add feedback control
                        t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                        u_fdb_k = u(:,k) - K*(Sig_init(1:obj.nx,j,k) - x(:,k));
                        
                        if k > 1
                            x_previous = xdn_previous_full(:,j);
                            df_previous = df_previous_full(:,:,j);
                        end
                        
                        % estimate whether current state is close to contact
                        [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig_init(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                        phi_bottom = phi_current(2:2:end);
                        active_threshold = 0.5;
                        if 0 %any(phi_bottom<active_threshold)
                            for kk = 1:size_terrain_sample
                                if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                    obj.plant.friction_coeff = w_mu(kk);
                                elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                    obj.plant.terrain_index = kk;
                                end
                                [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                            end
                            xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                            Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                            Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                        else
                            [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig_init(1:nx,j,k),u_fdb_k);
                            
                            % %numerical diff
                            % dt = diag(sqrt(eps(t)));
                            % dx = diag(max(sqrt(eps(Sig(1:obj.nx,j,k))), 1e-7*ones(obj.nx,1)));
                            % du = diag(max(sqrt(eps(u_fdb_k)), 1e-7*ones(nu,1)));
                            %
                            % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                            %
                            % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                            % for m = 1:N_finite_diff_x
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                            %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                            % end
                            %
                            % N_finite_diff_u = length(u_fdb_k);
                            % for m = 1:N_finite_diff_u
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                            %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                            % end
                            %
                            % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                            %     keyboard
                            % end
                            
                            xdn = xdn_analytical;
                            df(:,:,j) = df_analytical(:,:,j);
                            
                            Sig(1:nx,j,k+1) = xdn(1:nx,j);
                            dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                            dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                            dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                        end
                    end
                    xdn_previous_full = xdn;
                    df_previous_full = df;
                    
                    %calculate mean and variance w.r.t. [x_k] from sigma points
                    x_mean(:,k+1) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    end
                    
                    % % shift sigma point
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                    %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                    %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                    %         eta_final(k+1) = eta_optimal;
                    %     end
                    % end
                    %
                    % % scale Sigma point deviation by eta
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                    % end
                    %
                    % % recalculate mean and variance w.r.t. [x_k] from sigma points
                    % x_mean(:,k+1) = zeros(obj.nx,1);
                    % for j = 1:n_sampling_point
                    %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    % end
                    %
                    % % check that the mean deviation term is cancelled out
                    % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                    %     disp('shifting scheme is not correct')
                    %     keyboard
                    % end
                    
                    Px(:,:,k+1) = zeros(obj.nx);
                    for jj = 1:n_sampling_point
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    end
                    %c = c + kappa*trace(Px(:,:,k+1));
                    
                    % for jj = 1:n_sampling_point
                    % V_comp(:,:,k+1) = (Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    % c = c + kappa*trace(w*V_comp(:,:,k+1));
                    % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                    
                    % accumulate returned cost
                    %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                    c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*(x(:,k+1)-x_mean(:,k+1));%ML mean deviation version
                    
                    % derivative of variance matrix
                    % gradient of Tr(V) w.r.t state vector x
                    dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
                    dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                    
                    tic
                    % gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(obj.nx);
                    % gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                    dSig_i_kplus1_dx_set = zeros(obj.nx, obj.nx, n_sampling_point);
                    dSig_i_kplus1_du_set = zeros(obj.nx, nu, n_sampling_point);
                    
                    for i=1:n_sampling_point
                        if i == 1
                            for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                % gradient of Tr(V_{k+1}) w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(obj.nx);
                                dSig_m_kplus1_du = zeros(obj.nx,1);
                                
                                dSig_m_kplus1_dx = dfdx(:,:,m,k+1) + dfdSig(:,:,m,k+1);
                                dSig_m_kplus1_du = dfdu(:,:,m,k+1);% [double check that du is not affected]
                                
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run 2*(obj.nx+nw) times in total to obtain
                        % gradient w.r.t sigma points
                        dSig_i_kplus1_dx = zeros(obj.nx);
                        dSig_i_kplus1_du = zeros(obj.nx,nu);
                        
                        dSig_i_kplus1_dx = dfdx(:,:,i,k+1) + dfdSig(:,:,i,k+1);
                        dSig_i_kplus1_du = dfdu(:,:,i,k+1);
                        
                        dSig_i_kplus1_dx_set(:,:,i) = dSig_i_kplus1_dx;
                        dSig_i_kplus1_du_set(:,:,i) = dSig_i_kplus1_du;
                    end
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    dmeanRdx(k,:,k+1) = dmeanRdx(k,:,k+1) + ((x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    dmeanRdu(k,:,k+1) = dmeanRdu(k,:,k+1) + ((x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    
                    dCovdx = zeros(nx,nx,nx);
                    dCovdu = zeros(nx,nx,nu);
                    
                    %                 for pp=1:nx
                    %                     for i=1:n_sampling_point
                    %                         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                    %                         for nn=1:nx
                    %                             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                    %                             if nn == 1
                    %                                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                    %                                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                    %                             elseif nn == nx
                    %                                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                    %                             else
                    %                                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                    %                                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                    %                             end
                    %                             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_set(nn,pp,i) - w_avg*dSig_m_kplus1_dx_sum(nn,pp));
                    %                             if pp <= nu
                    %                                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_set(nn,pp,i) - w_avg*dSig_m_kplus1_du_sum(nn,pp));
                    %                             end
                    %                         end
                    %                     end
                    %                     dmeanRdx(k,pp,k+1) = dmeanRdx(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                    %                     dmeanRdx(k,pp,k+1) = dmeanRdx(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                    %                     if pp <= nu
                    %                         dmeanRdu(k,pp,k+1) = dmeanRdu(k,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                    %                         dmeanRdu(k,pp,k+1) = dmeanRdu(k,pp,k+1) + 1/2*trace(pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                    %                     end
                    %                 end
                end
                
                dc = [];
                % cost gradient w.r.t x at first time step is zero
                for jj=1:obj.N % index for x_k
                    if (jj == 1)
                        dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                    else
                        %dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                        dmeanR_sum_dx_k = ((x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                    end
                    
                    if jj < obj.N
                        dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jj+1);
                    end
                    %dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                    dc = [dc, dmeanR_sum_dx_k];
                end
                
                % cost gradient w.r.t u at first time step is zero, since
                for jj=1:obj.N % index for u_k
                    dmeanR_sum_du_k = zeros(1, nu);
                    if jj < obj.N
                        dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jj+1);
                    end
                    %dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                    dc = [dc, dmeanR_sum_du_k];
                end
                
                % scale this robust cost
                c = obj.options.contact_robust_cost_coeff*c;
                dc = obj.options.contact_robust_cost_coeff*dc;
                
            end
        end
        
        function [c,dc] = robustVariancecost(obj, x_full, Fext_full)
            global timestep_updated
            global time_step
            global x_previous
            global df_previous
            global terrain_index
            global uncertain_mu
            manip = obj.plant.getManipulator();
            
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.nFext;%obj.plant.getNumInputs;
            obj.nu = nu;
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px(:,:,1) = obj.cached_Px(:,:,1);
            
            if strcmp(manip.uncertainty_source,'friction_coeff')
                w_mu = manip.uncertain_mu_set;
                w_noise = [w_mu];
                Pw = diag([0.01]);
            elseif strcmp(manip.uncertainty_source,'terrain_height')
                w_phi = manip.uncertain_position_set;
                w_noise = [w_phi];
                Pw = diag([0.016]);
            elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                w_mu = manip.uncertain_mu_set;
                w_phi = manip.uncertain_position_set;
                w_noise = [w_mu;w_phi];
                Pw = diag([0.01, 0.016]);
            elseif isempty(manip.uncertainty_source)
                Pw = [];
                w_noise = [];
            elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                save -ascii friction_coeff_noise.dat w_mu
                %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                %w_phi = [x;y];%height noise
                %save -ascii initial_position_noise.dat w_phi
                w_noise = [w_mu];
            end
            
            % disturbance variance
            % currently only consider terrain height and/or friction coefficient
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            nw = size(Pw,1);
            n_sampling_point = 2*(obj.nx+nw);
            w_avg = 1/n_sampling_point;
            
            K = obj.options.K;
            eta_final = ones(obj.N,1);
            
            %initialize c and dc
            kappa = obj.options.kappa;
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            c = 0;
            c = kappa*trace(obj.cached_Px(:,:,1));
            %c_quadratic = 0;
            %c_variance = 0;
            dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
            
            % initialize gradient of Tr(V) w.r.t state vector x
            dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
            dTrVdu(:,:,1) = zeros(obj.N-1,nu);
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
            end
            
            plant_update = @objPlantUpdate;
            %obj.plant.time_step = time_step;
            df_previous_full = [];
            xdn_previous_full = [];
            
            for k = 1:obj.N-1%[Ye: double check the index]
                %Propagate sigma points through nonlinear dynamics
                if k == 1
                    [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                    if d
                        diverge = k;
                        return;
                    end
                    %S = zeros(length(Px(:,:,k)));
                    
                    S = scale*S;
                    Sig(:,:,k) = [S -S];
                    if isempty(w_noise)
                        for j = 1:n_sampling_point
                            Sig(:,j,k) = Sig(:,j,k) + [x(:,k)];
                        end
                    else
                        for j = 1:n_sampling_point
                            Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                        end
                    end
                    x_mean(:,k) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:obj.nx,j,k);
                    end
                    c = c + norm(x(:,k)-x_mean(:,k))^2;
                    
                    % % debugging
                    % c_quadratic(k) = norm(x(:,k)-x_mean(:,k))^2;
                    % c_quadratic_x(k) = norm(x(1,k)-x_mean(1,k))^2;
                    % c_quadratic_xd(k) = norm(x(7,k)-x_mean(7,k))^2;
                    % c_quadratic_z(k) = norm(x(3,k)-x_mean(3,k))^2;
                    % c_quadratic_zd(k) = norm(x(9,k)-x_mean(9,k))^2;
                    % for j = 1:n_sampling_point
                    %     V_comp = (Sig(1:obj.nx,j,k)-x_mean(:,k))*(Sig(1:obj.nx,j,k)-x_mean(:,k))';
                    %     c_variance(j,k) = kappa*trace(w*V_comp);
                    %
                    %     % debugging
                    %     V_comp_x = (Sig(1,j,k)-x_mean(1,k))*(Sig(1,j,k)-x_mean(1,k))';
                    %     c_variance_x(j,k) = kappa*trace(w*V_comp_x);
                    %     V_comp_xd = (Sig(7,j,k)-x_mean(7,k))*(Sig(7,j,k)-x_mean(7,k))';
                    %     c_variance_xd(j,k) = kappa*trace(w*V_comp_xd);
                    %
                    %     V_comp_z = (Sig(3,j,k)-x_mean(3,k))*(Sig(3,j,k)-x_mean(3,k))';
                    %     c_variance_z(j,k) = kappa*trace(w*V_comp_z);
                    %     V_comp_zd = (Sig(9,j,k)-x_mean(9,k))*(Sig(9,j,k)-x_mean(9,k))';
                    %     c_variance_zd(j,k) = kappa*trace(w*V_comp_zd);
                    % end
                end
                
                %Propagate sigma points through nonlinear dynamics
                for j = 1:n_sampling_point
                    noise_index = j;
                    
                    % a hacky way to implement the control input
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = w_mu(j);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        terrain_index = j;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        uncertain_mu = w_mu(j);
                        terrain_index = j;
                    end
                    
                    % add feedback control
                    t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                    
                    if k > 1
                        x_previous = xdn_previous_full(:,j);
                        df_previous = df_previous_full(:,:,j);
                    end
                    
                    % estimate whether current state is close to contact
                    [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                    phi_bottom = phi_current(2:2:end);
                    active_threshold = 0.5;
                    if 0 %any(phi_bottom<active_threshold)
                        for kk = 1:size_terrain_sample
                            if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                obj.plant.friction_coeff = w_mu(kk);
                            elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                obj.plant.terrain_index = kk;
                            end
                            [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                        end
                        xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                        Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                        Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                        
                    else
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig(1:nx,j,k),u_fdb_k);
                        
                        % %numerical diff
                        % dt = diag(sqrt(eps(t)));
                        % dx = diag(max(sqrt(eps(Sig(1:obj.nx,j,k))), 1e-3*ones(obj.nx,1)));
                        % du = diag(max(sqrt(eps(u_fdb_k)), 1e-3*ones(nu,1)));
                        %
                        % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                        %
                        % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                        % for m = 1:N_finite_diff_x
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                        %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                        % end
                        %
                        % N_finite_diff_u = length(u_fdb_k);
                        % for m = 1:N_finite_diff_u
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                        %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                        % end
                        %
                        % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                        %     keyboard
                        % end
                        
                        xdn = xdn_analytical;
                        df(:,:,j) = df_analytical(:,:,j);
                        
                        Sig(1:nx,j,k+1) = xdn(1:nx,j);
                        dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                        dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                    end
                end
                xdn_previous_full = xdn;
                df_previous_full = df;
                
                %calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(obj.nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                end
                
                % % shift sigma point
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                %         eta_final(k+1) = eta_optimal;
                %     end
                % end
                %
                % % scale Sigma point deviation by eta
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                % end
                %
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(obj.nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                % end
                %
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                
                Px(:,:,k+1) = zeros(obj.nx);
                for jj = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                end
                c = c + kappa*trace(Px(:,:,k+1));
                
                % for jj = 1:n_sampling_point
                % V_comp(:,:,k+1) = (Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                % c = c + kappa*trace(w*V_comp(:,:,k+1));
                % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                
                % % debugging
                % c_variance(j,k+1) = kappa*trace(w*V_comp);
                %
                % V_comp_x = (Sig(1,j,k+1)-x_mean(1,k+1))*(Sig(1,j,k+1)-x_mean(1,k+1))';
                % c_variance_x(j,k+1) = kappa*trace(w*V_comp_x);
                % V_comp_xd = (Sig(7,j,k+1)-x_mean(7,k+1))*(Sig(7,j,k+1)-x_mean(7,k+1))';
                % c_variance_xd(j,k+1) = kappa*trace(w*V_comp_xd);
                %
                % V_comp_z = (Sig(3,j,k+1)-x_mean(3,k+1))*(Sig(3,j,k+1)-x_mean(3,k+1))';
                % c_variance_z(j,k+1) = kappa*trace(w*V_comp_z);
                % V_comp_zd = (Sig(9,j,k+1)-x_mean(9,k+1))*(Sig(9,j,k+1)-x_mean(9,k+1))';
                % c_variance_zd(j,k+1) = kappa*trace(w*V_comp_zd);
                % end
                
                % accumulate returned cost
                c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                % % debugging
                % c_quadratic(k+1) = norm(x(:,k+1)-x_mean(:,k+1))^2;
                % c_quadratic_x(k+1) = norm(x(1,k+1)-x_mean(1,k+1))^2;
                % c_quadratic_xd(k+1) = norm(x(7,k+1)-x_mean(7,k+1))^2;
                % c_quadratic_z(k+1) = norm(x(3,k+1)-x_mean(3,k+1))^2;
                % c_quadratic_zd(k+1) = norm(x(7,k+1)-x_mean(9,k+1))^2;
                
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
                    
                    tic
                    % gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(obj.nx);
                    % gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                    dSig_i_kplus1_dx_resample = zeros(obj.nx,obj.nx,n_sampling_point);
                    dSig_i_kplus1_du_resample = zeros(obj.nx,nu,n_sampling_point);
                    
                    for i=1:n_sampling_point
                        if i == 1
                            for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                % gradient of Tr(V_{k+1}) w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(obj.nx);
                                dSig_m_kplus1_du = zeros(obj.nx,1);
                                
                                chain_rule_indx = k-j;
                                if j ~= 1
                                    %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                    dSig_m_kplus1_dx = dfdx(:,:,m,2) + dfdSig(:,:,m,2);
                                end
                                dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]

                                while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                    chain_rule_indx = chain_rule_indx - 1;
                                end

                                % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
                                %     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
                                %     chain_rule_indx = chain_rule_indx - 1;
                                % end
                                %
                                % if (chain_rule_indx > 0)
                                %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                                %         chain_rule_indx = chain_rule_indx - 1;
                                %     end
                                % else
                                %     chain_rule_grad = eye(12);
                                % end
                                %
                                % dSig_m_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_dx;
                                % dSig_m_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_du;
                                
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run 2*(obj.nx+nw) times in total to obtain
                        % gradient w.r.t sigma points
                        dSig_i_kplus1_dx = zeros(obj.nx);
                        dSig_i_kplus1_du = zeros(obj.nx,nu);
                        chain_rule_indx = k-j;
                        if j ~= 1
                            %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                        else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            dSig_i_kplus1_dx = dfdx(:,:,i,2) + dfdSig(:,:,i,2)*eye(obj.nx);
                        end
                        dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                            dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                        % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        %     dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_dx;
                        %     dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_du;
                        %     chain_rule_indx = chain_rule_indx - 1;
                        % end
                        %
                        % if (chain_rule_indx > 0)
                        %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                        %         chain_rule_indx = chain_rule_indx - 1;
                        %     end
                        % else
                        %     chain_rule_grad = eye(obj.nx);
                        % end
                        %
                        % dSig_i_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_dx;
                        % dSig_i_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_du;
                        
                          dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                          dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
                          % % new sigma point due to resampling mechanism
                          % dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                          % dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                    end
                    
                    % %% resample scheme
                    % dSig_kplus1_dx_sum_resample = zeros(obj.nx,obj.nx);
                    % dSig_kplus1_du_sum_resample = zeros(obj.nx,nu);
                    %
                    % for i =1:n_sampling_point
                    %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                    %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                    % end
                    %
                    % for i =1:n_sampling_point
                    %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                    %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                    % end
                    % %% resample scheme
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                    %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    
                    % dCovdx = zeros(nx,nx,nx);
                    % dCovdu = zeros(nx,nx,nu);
                    % toc
                    %
                    % tic
                    % for pp=1:nx
                    %     for i=1:n_sampling_point
                    %         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                    %         for nn=1:nx
                    %             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                    %             if nn == 1
                    %                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                    %                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                    %             elseif nn == nx
                    %                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                    %             else
                    %                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                    %                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                    %             end
                    %             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_resample(nn,pp,i) - w_avg*dSig_kplus1_dx_sum_resample(nn,pp));
                    %             if pp <= nu
                    %                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_resample(nn,pp,i) - w_avg*dSig_kplus1_du_sum_resample(nn,pp));
                    %             end
                    %         end
                    %     end
                    %     dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                    %     if pp <= nu
                    %         dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                    %     end
                    % end
                    % toc
                end
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
                dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                %dc = [dc, dmeanR_sum_dx_k];
                %dc = [dc, kappa*dTrV_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            for jj=1:obj.N % index for u_k
                dTrV_sum_du_k = zeros(1, nu);
                dmeanR_sum_du_k = zeros(1, nu);
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                end
                dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                %dc = [dc, dmeanR_sum_du_k];
                %dc = [dc, kappa*dTrV_sum_du_k];
            end
            
            % scale this robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
            fprintf('robust cost function: %4.8f\n',c);
            tElapsed = toc(tStart);
            
            % figure(7),hold on;plot(c_quadratic_x,'b-')
            % figure(8),hold on;plot(c_quadratic_xd,'b-')
            % figure(9),hold on;plot(c_variance_x(1,:),'b-')
            % figure(10),hold on;plot(c_variance_xd(1,:),'b-')
            %
            % figure(11),hold on;plot(c_quadratic_z,'b-')
            % figure(12),hold on;plot(c_quadratic_zd,'b-')
            % figure(13),hold on;plot(c_variance_z(1,:),'b-')
            % figure(14),hold on;plot(c_variance_zd(1,:),'b-')
            %
            % figure(15)
            % clf
            % Sig_permute = permute(Sig,[1,3,2]);
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(1,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(1,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,7])
            %
            % figure(16)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(7,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(7,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,11])
            %
            % figure(17)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(3,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(3,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,4])
            %
            % figure(18)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(9,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(9,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([-7,0])
            %
            % figure(19)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(5,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(5,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([-1,1])
            %
            % obj.cached_Px = Px;
            % fprintf('robust cost function: %4.8f\n',c);
            
            % % check gradient
            % disp('check gradient')
            % c_numeric = c;
            % dc_numeric = dc;
            %
            X0 = [x_full; Fext_full];
            % %X0 = X0 + randn(size(X0))*0.1;
            %
            % fun = @(X0) robustVariancecost_check(obj, X0);
            % DerivCheck(fun, X0)
            % disp('finish')
            c_numeric = c;
            dc_numeric = dc;
             
            % [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_check(obj,X0),X0,struct('grad_method','numerical'));
            %
            % [c_numeric,dc_numeric] = robustVariancecost_check(obj, X0);
            %
            % valuecheck(dc,dc_numeric,1e-5);
            % valuecheck(c,c_numeric,1e-5);
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %
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
                dd =  dot(rr, JJ)-f2+f1;
                fprintf('difference between numerical and analytical: %4.15f\n',dd);
            end
            
            %             function [c,dc] = robustVariancecost_check(obj, X0)
            %                 x_full = X0(1:obj.nx*obj.N);
            %                 Fext_full = X0(obj.nx*obj.N+1:end);
            %
            %                 x = reshape(x_full, obj.nx, obj.N);
            %                 u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            %                 nq = obj.plant.getNumPositions;
            %                 nv = obj.plant.getNumVelocities;
            %                 nu = obj.nFext;%obj.plant.getNumInputs;
            
            function [c,dc] = robustVariancecost_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                
                manip = obj.plant.getManipulator();
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nx = nq+nv;
                nu = obj.nFext;%obj.plant.getNumInputs;
                obj.nu = nu;
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px(:,:,1) = obj.cached_Px(:,:,1);
                
                if strcmp(manip.uncertainty_source,'friction_coeff')
                    w_mu = manip.uncertain_mu_set;
                    w_noise = [w_mu];
                    Pw = diag([0.01]);
                elseif strcmp(manip.uncertainty_source,'terrain_height')
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_phi];
                    Pw = diag([0.016]);
                elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                    w_mu = manip.uncertain_mu_set;
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_mu;w_phi];
                    Pw = diag([0.01, 0.016]);
                elseif isempty(manip.uncertainty_source)
                    Pw = [];
                    w_noise = [];
                elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                    %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                    w_noise = [w_mu];
                end
                
                % disturbance variance
                % currently only consider terrain height and/or friction coefficient
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
                nw = size(Pw,1);
                n_sampling_point = 2*(obj.nx+nw);
                w_avg = 1/n_sampling_point;
                
                K = obj.options.K;
                eta_final = ones(obj.N,1);
                
                %initialize c and dc
                kappa = obj.options.kappa;
                x_mean = zeros(obj.nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
                c = 0;
                c = kappa*trace(obj.cached_Px(:,:,1));
                %c_quadratic = 0;
                %c_variance = 0;
                dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
                
                % initialize gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
                dTrVdu(:,:,1) = zeros(obj.N-1,nu);
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                    [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
                end
                
                plant_update = @objPlantUpdate;
                %obj.plant.time_step = time_step;
                df_previous_full = [];
                xdn_previous_full = [];
                
                for k = 1:obj.N-1%[Ye: double check the index]
                    %Propagate sigma points through nonlinear dynamics
                    if k == 1
                        [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                        if d
                            diverge = k;
                            return;
                        end
                        %S = zeros(length(Px(:,:,k)));
                        
                        S = scale*S;
                        Sig(:,:,k) = [S -S];
                        if isempty(w_noise)
                            for j = 1:n_sampling_point
                                Sig(:,j,k) = Sig(:,j,k) + [x(:,k)];
                            end
                        else
                            for j = 1:n_sampling_point
                                Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                            end
                        end
                        x_mean(:,k) = zeros(obj.nx,1);
                        for j = 1:n_sampling_point
                            x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:obj.nx,j,k);
                        end
                        c = c + norm(x(:,k)-x_mean(:,k))^2;
                        
                        % % debugging
                        % c_quadratic(k) = norm(x(:,k)-x_mean(:,k))^2;
                        % c_quadratic_x(k) = norm(x(1,k)-x_mean(1,k))^2;
                        % c_quadratic_xd(k) = norm(x(7,k)-x_mean(7,k))^2;
                        % c_quadratic_z(k) = norm(x(3,k)-x_mean(3,k))^2;
                        % c_quadratic_zd(k) = norm(x(9,k)-x_mean(9,k))^2;
                        % for j = 1:n_sampling_point
                        %     V_comp = (Sig(1:obj.nx,j,k)-x_mean(:,k))*(Sig(1:obj.nx,j,k)-x_mean(:,k))';
                        %     c_variance(j,k) = kappa*trace(w*V_comp);
                        %
                        %     % debugging
                        %     V_comp_x = (Sig(1,j,k)-x_mean(1,k))*(Sig(1,j,k)-x_mean(1,k))';
                        %     c_variance_x(j,k) = kappa*trace(w*V_comp_x);
                        %     V_comp_xd = (Sig(7,j,k)-x_mean(7,k))*(Sig(7,j,k)-x_mean(7,k))';
                        %     c_variance_xd(j,k) = kappa*trace(w*V_comp_xd);
                        %
                        %     V_comp_z = (Sig(3,j,k)-x_mean(3,k))*(Sig(3,j,k)-x_mean(3,k))';
                        %     c_variance_z(j,k) = kappa*trace(w*V_comp_z);
                        %     V_comp_zd = (Sig(9,j,k)-x_mean(9,k))*(Sig(9,j,k)-x_mean(9,k))';
                        %     c_variance_zd(j,k) = kappa*trace(w*V_comp_zd);
                        % end
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:n_sampling_point
                        noise_index = j;
                        
                        % a hacky way to implement the control input
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        
                        if strcmp(manip.uncertainty_source, 'friction_coeff')
                            uncertain_mu = w_mu(j);
                        elseif strcmp(manip.uncertainty_source, 'terrain_height')
                            terrain_index = j;
                        elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                            uncertain_mu = w_mu(j);
                            terrain_index = j;
                        end
                        
                        % add feedback control
                        t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                        u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                        
                        if k > 1
                            x_previous = xdn_previous_full(:,j);
                            df_previous = df_previous_full(:,:,j);
                        end
                        
                        % estimate whether current state is close to contact
                        [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                        phi_bottom = phi_current(2:2:end);
                        active_threshold = 0.5;
                        if 0 %any(phi_bottom<active_threshold)
                            for kk = 1:size_terrain_sample
                                if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                    obj.plant.friction_coeff = w_mu(kk);
                                elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                    obj.plant.terrain_index = kk;
                                end
                                [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                            end
                            xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                            Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                            Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                            
                        else
                            [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig(1:nx,j,k),u_fdb_k);
                            
                            % %numerical diff
                            % dt = diag(sqrt(eps(t)));
                            % dx = diag(max(sqrt(eps(Sig(1:obj.nx,j,k))), 1e-3*ones(obj.nx,1)));
                            % du = diag(max(sqrt(eps(u_fdb_k)), 1e-3*ones(nu,1)));
                            %
                            % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                            %
                            % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                            % for m = 1:N_finite_diff_x
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                            %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                            % end
                            %
                            % N_finite_diff_u = length(u_fdb_k);
                            % for m = 1:N_finite_diff_u
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                            %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                            % end
                            %
                            % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                            %     keyboard
                            % end
                            
                            xdn = xdn_analytical;
                            df(:,:,j) = df_analytical(:,:,j);
                            
                            Sig(1:nx,j,k+1) = xdn(1:nx,j);
                            dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                            dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                            dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                        end
                    end
                    xdn_previous_full = xdn;
                    df_previous_full = df;
                    
                    %calculate mean and variance w.r.t. [x_k] from sigma points
                    x_mean(:,k+1) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    end
                    
                    % % shift sigma point
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                    %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                    %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                    %         eta_final(k+1) = eta_optimal;
                    %     end
                    % end
                    %
                    % % scale Sigma point deviation by eta
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                    % end
                    %
                    % % recalculate mean and variance w.r.t. [x_k] from sigma points
                    % x_mean(:,k+1) = zeros(obj.nx,1);
                    % for j = 1:n_sampling_point
                    %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    % end
                    %
                    % % check that the mean deviation term is cancelled out
                    % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                    %     disp('shifting scheme is not correct')
                    %     keyboard
                    % end
                    
                    Px(:,:,k+1) = zeros(obj.nx);
                    for jj = 1:n_sampling_point
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    end
                    c = c + kappa*trace(Px(:,:,k+1));
                    
                    % for jj = 1:n_sampling_point
                    % V_comp(:,:,k+1) = (Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    % c = c + kappa*trace(w*V_comp(:,:,k+1));
                    % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                    
                    % % debugging
                    % c_variance(j,k+1) = kappa*trace(w*V_comp);
                    %
                    % V_comp_x = (Sig(1,j,k+1)-x_mean(1,k+1))*(Sig(1,j,k+1)-x_mean(1,k+1))';
                    % c_variance_x(j,k+1) = kappa*trace(w*V_comp_x);
                    % V_comp_xd = (Sig(7,j,k+1)-x_mean(7,k+1))*(Sig(7,j,k+1)-x_mean(7,k+1))';
                    % c_variance_xd(j,k+1) = kappa*trace(w*V_comp_xd);
                    %
                    % V_comp_z = (Sig(3,j,k+1)-x_mean(3,k+1))*(Sig(3,j,k+1)-x_mean(3,k+1))';
                    % c_variance_z(j,k+1) = kappa*trace(w*V_comp_z);
                    % V_comp_zd = (Sig(9,j,k+1)-x_mean(9,k+1))*(Sig(9,j,k+1)-x_mean(9,k+1))';
                    % c_variance_zd(j,k+1) = kappa*trace(w*V_comp_zd);
                    % end
                    
                    % accumulate returned cost
                    c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                    % % debugging
                    % c_quadratic(k+1) = norm(x(:,k+1)-x_mean(:,k+1))^2;
                    % c_quadratic_x(k+1) = norm(x(1,k+1)-x_mean(1,k+1))^2;
                    % c_quadratic_xd(k+1) = norm(x(7,k+1)-x_mean(7,k+1))^2;
                    % c_quadratic_z(k+1) = norm(x(3,k+1)-x_mean(3,k+1))^2;
                    % c_quadratic_zd(k+1) = norm(x(7,k+1)-x_mean(9,k+1))^2;
                    
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
                        
                        tic
                        % gradient w.r.t state x
                        dSig_m_kplus1_dx_sum = zeros(obj.nx);
                        % gradient w.r.t control u
                        dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                        dSig_i_kplus1_dx_resample = zeros(obj.nx,obj.nx,n_sampling_point);
                        dSig_i_kplus1_du_resample = zeros(obj.nx,nu,n_sampling_point);
                        
                        for i=1:n_sampling_point
                            if i == 1
                                for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                    % gradient of Tr(V_{k+1}) w.r.t control x and u
                                    dSig_m_kplus1_dx = zeros(obj.nx);
                                    dSig_m_kplus1_du = zeros(obj.nx,1);
                                    
                                    chain_rule_indx = k-j;
                                    if j ~= 1
                                        %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                    else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                        dSig_m_kplus1_dx = dfdx(:,:,m,2) + dfdSig(:,:,m,2);
                                    end
                                    dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                    
                                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                        dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                        dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                        chain_rule_indx = chain_rule_indx - 1;
                                    end
                                    
                                    % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    %     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
                                    %     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
                                    %     chain_rule_indx = chain_rule_indx - 1;
                                    % end
                                    %
                                    % if (chain_rule_indx > 0)
                                    %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                                    %         chain_rule_indx = chain_rule_indx - 1;
                                    %     end
                                    % else
                                    %     chain_rule_grad = eye(12);
                                    % end
                                    %
                                    % dSig_m_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_dx;
                                    % dSig_m_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_du;
                                    
                                    dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                                end
                            end
                            
                            % run 2*(obj.nx+nw) times in total to obtain
                            % gradient w.r.t sigma points
                            dSig_i_kplus1_dx = zeros(obj.nx);
                            dSig_i_kplus1_du = zeros(obj.nx,nu);
                            chain_rule_indx = k-j;
                            if j ~= 1
                                %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                            else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                dSig_i_kplus1_dx = dfdx(:,:,i,2) + dfdSig(:,:,i,2)*eye(obj.nx);
                            end
                            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            
                            % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            %     dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_dx;
                            %     dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_du;
                            %     chain_rule_indx = chain_rule_indx - 1;
                            % end
                            %
                            % if (chain_rule_indx > 0)
                            %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                            %         chain_rule_indx = chain_rule_indx - 1;
                            %     end
                            % else
                            %     chain_rule_grad = eye(obj.nx);
                            % end
                            %
                            % dSig_i_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_dx;
                            % dSig_i_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_du;
                            
                            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
                            % % new sigma point due to resampling mechanism
                            % dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            % dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        end
                        
                        % %% resample scheme
                        % dSig_kplus1_dx_sum_resample = zeros(obj.nx,obj.nx);
                        % dSig_kplus1_du_sum_resample = zeros(obj.nx,nu);
                        %
                        % for i =1:n_sampling_point
                        %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                        %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                        % end
                        %
                        % for i =1:n_sampling_point
                        %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                        %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                        % end
                        % %% resample scheme
                        
                        % gradient of mean residual w.r.t state x and control u, assume norm 2
                        dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                        dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                        %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                        %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                        
                        % dCovdx = zeros(nx,nx,nx);
                        % dCovdu = zeros(nx,nx,nu);
                        % toc
                        %
                        % tic
                        % for pp=1:nx
                        %     for i=1:n_sampling_point
                        %         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                        %         for nn=1:nx
                        %             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                        %             if nn == 1
                        %                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                        %                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                        %             elseif nn == nx
                        %                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                        %                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                        %             else
                        %                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                        %                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                        %                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                        %             end
                        %             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_resample(nn,pp,i) - w_avg*dSig_kplus1_dx_sum_resample(nn,pp));
                        %             if pp <= nu
                        %                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_resample(nn,pp,i) - w_avg*dSig_kplus1_du_sum_resample(nn,pp));
                        %             end
                        %         end
                        %     end
                        %     dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                        %     if pp <= nu
                        %         dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                        %     end
                        % end
                        % toc
                    end
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
                    dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                    %dc = [dc, dmeanR_sum_dx_k];
                    %dc = [dc, kappa*dTrV_sum_dx_k];
                end
                
                % cost gradient w.r.t u at first time step is zero, since
                for jj=1:obj.N % index for u_k
                    dTrV_sum_du_k = zeros(1, nu);
                    dmeanR_sum_du_k = zeros(1, nu);
                    for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                        dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                    end
                    dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                    %dc = [dc, dmeanR_sum_du_k];
                    %dc = [dc, kappa*dTrV_sum_du_k];
                end
                
                % scale this robust cost
                c = obj.options.contact_robust_cost_coeff*c;
                dc = obj.options.contact_robust_cost_coeff*dc;
            end
            
            function eta_optimal = etaEstimation(xdn,x_nominal)
                % hard coding limit for falling brick position and orientation
                % q_lb and q_ub(9:11) to be defined
                
                eta_optimal = 1;% non-scaling
                for i = 1:length(q_lb)
                    eta = 0;
                    if xdn(i) < q_lb(i)
                        eta = abs(q_lb(i)/xdn(i));
                    elseif xdn(i) > q_ub(i)
                        eta = abs(q_ub(i)/xdn(i));
                    end
                    if eta ~= 0 && eta < eta_optimal
                        eta_optimal = eta;
                    end
                end
                
                % limit on velocity
                v_lb = [-2.5*ones(7,1);-0.2;-2*ones(3,1);-2.5*ones(3,1)];
                v_ub = [2.5*ones(7,1);0.2;2*ones(3,1);2.5*ones(3,1)];
                for i = 1:length(v_lb)
                    eta = 0;
                    if xdn(nq+i) < v_lb(i)
                        eta = abs(v_lb(i)/xdn(nq+i));
                    elseif xdn(nq+i) > v_ub(i)
                        eta = abs(v_ub(i)/xdn(nq+i));
                    end
                    if eta ~= 0 && eta < eta_optimal
                        eta_optimal = eta;
                    end
                end
                eta_optimal = eta_optimal/5;% scale down
            end
        end
        
        function [c,dc] = robustVariancecost_scaled(obj, x_full, Fext_full)
            global timestep_updated
            global x_previous
            global df_previous
            global terrain_index
            global uncertain_mu
            manip = obj.plant.getManipulator();
             
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.nFext;%obj.plant.getNumInputs;
            obj.nu = nu;
            
            % sigma points
            Px = zeros(nx,nx,obj.N);
            Px(:,:,1) = obj.cached_Px(:,:,1);
            Px_init = Px(:,:,1);
            
            if strcmp(manip.uncertainty_source,'friction_coeff')
                w_mu = manip.uncertain_mu_set;
                w_noise = [w_mu];
                Pw = diag([0.01]);
            elseif strcmp(manip.uncertainty_source,'terrain_height')
                w_phi = manip.uncertain_position_set;
                w_noise = [w_phi];
                Pw = diag([0.016]);
            elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                w_mu = manip.uncertain_mu_set;
                w_phi = manip.uncertain_position_set;
                w_noise = [w_mu;w_phi];
                Pw = diag([0.01, 0.016]);
            elseif isempty(manip.uncertainty_source)
                Pw = [];
                w_noise = [];
            elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                save -ascii friction_coeff_noise.dat w_mu
                %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                %w_phi = [x;y];%height noise
                %save -ascii initial_position_noise.dat w_phi
                w_noise = [w_mu];
            end
            
            % disturbance variance
            % currently only consider terrain height and/or friction coefficient
            
            nw = size(Pw,1);
            sampling_method = 1;%option 1: unscented transform, option 2: random sampling with a smaller number
            if sampling_method == 1
                n_sampling_point = 2*(nx+nw);
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
            elseif sampling_method == 2
                n_sampling_point = 10;
                w_state = load('state_noise.dat');%std:0.001
                %w_state = 0.001*randn(28,62);
                %save -ascii state_noise_small.dat w_state
                w = 1;
            end
            w_avg = 1/n_sampling_point;
            K = obj.options.K;
            alpha = obj.options.alpha;
            
            %initialize c and dc
            kappa = obj.options.kappa;
            x_mean = zeros(nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            c = 0;
            c_quadratic = 0;
            c_covariance = 0;
            dc = zeros(1, 1+obj.N*(nx+nu));
            
            % initialize gradient of Tr(V) w.r.t state vector x
            dTrVdx(:,:,1) = zeros(obj.N-1,nx);
            dTrVdu(:,:,1) = zeros(obj.N-1,nu);
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(timestep_updated,Sig,u_fdb_k)
                [xdn,df] = obj.plant.update(timestep_updated,Sig,u_fdb_k);
            end
            
            plant_update = @objPlantUpdate;
            obj.plant.time_step = timestep_updated;
            df_previous_full = [];
            xdn_previous_full = [];
            noise_sample_type = 1;
            
            for k = 1:obj.N-1%[Ye: double check the index]
                if sampling_method == 1
                    %Propagate sigma points through nonlinear dynamics
                    if k == 1
                        [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
                        
                        if d
                            diverge = k;
                            return;
                        end
                        S = scale*S;
                        Sig_init(:,:,k) = [S -S];
                        
                        %state dimensionality of sigma point depends on the number of uncertainties.
                        %Our case w_noise has 0, 1, or 2 dimensions
                        if isempty(w_noise)
                            for j = 1:n_sampling_point
                                Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
                            end
                        else
                            for j = 1:n_sampling_point
                                Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
                            end
                        end
                        Sig(:,:,k) = Sig_init(:,:,k);
                        c = c + kappa*trace(Px_init);
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
                        c_quadratic = c_quadratic + 1/2*norm(x(:,k)-x_mean(:,k))^2;
                        
                        Px_init = zeros(nx);
                        for j = 1:n_sampling_point
                            Px_init = Px_init + w*(Sig(1:nx,j,k)-x_mean(:,k))*(Sig(1:nx,j,k)-x_mean(:,k))';
                        end
                        c = c + trace(Px_init);
                        c_covariance = c_covariance + trace(Px_init);
                    else
                        for j = 1:n_sampling_point
                            Sig_init(:,j,k) = (1-alpha)*Sig(:,j,k) + alpha*x(:,k);
                        end
                    end
                end
                
                %Propagate sigma points through nonlinear dynamics
                for j = 1:n_sampling_point                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = w_mu(j);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        terrain_index = j;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        uncertain_mu = w_mu(j);
                        terrain_index = j;
                    end
                    
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:nx/2,j,k),Sig_init(nx/2+1:nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    
                    % add feedback control
                    t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                    
                    if k > 1
                        x_previous = xdn_previous_full(:,j);
                        df_previous = df_previous_full(:,:,j);
                    end
                    
                    if noise_sample_type == 1
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                    elseif noise_sample_type == 2
                        xdn_analytical(:,j) = zeros(nx,1);
                        df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                        for kk=1:length(w_noise)
                            [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                            xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                            df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                        end
                        xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                        df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                    end
                    % %numerical diff
                    % dt = diag(max(sqrt(eps(timestep_updated)), 1e-7));
                    % dx = diag(max(sqrt(eps(Sig_init(1:nx,j,k))), 1e-7*ones(nx,1)));
                    % du = diag(max(sqrt(eps(u_fdb_k)), 1e-7*ones(nu,1)));
                    %
                    % [xdnp,~] = feval(plant_update,timestep_updated+dt,Sig_init(1:nx,j,k),u_fdb_k);
                    % [xdnm,~] = feval(plant_update,timestep_updated-dt,Sig_init(1:nx,j,k),u_fdb_k);
                    % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                    %
                    % N_finite_diff_x = length(Sig_init(1:nx,j,k));
                    % for m = 1:N_finite_diff_x
                    %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k)+dx(:,m),u_fdb_k);
                    %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k)-dx(:,m),u_fdb_k);
                    %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                    % end
                    %
                    % N_finite_diff_u = length(u_fdb_k);
                    % for m = 1:N_finite_diff_u
                    %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k+du(:,m));
                    %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k-du(:,m));
                    %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                    % end
                    
                    % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                    %     keyboard
                    % end
                    
                    xdn(:,j) = xdn_analytical(:,j);
                    df(:,:,j) = df_analytical(:,:,j);
                    
                    Sig(1:nx,j,k+1) = xdn(1:nx,j);
                    dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                    dfdSig(:,:,j,k+1) = (1-alpha)*df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*(1-alpha)*K;
                    dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*(1-alpha)*K + alpha*df(:,2:nx+1,j);
                end
                xdn_previous_full = xdn;
                df_previous_full = df;
                
                %calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                end
                
                % % shift sigma point
                % for j = 1:n_sampling_point
                %     Sig(1:nx,j,k+1) = Sig(1:nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                %     eta_optimal = 1;%etaEstimation(Sig(1:nx,j,k+1),x(:,k+1));
                %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                %         eta_final(k+1) = eta_optimal;
                %     end
                % end
                %
                % % scale Sigma point deviation by eta
                % for j = 1:n_sampling_point
                %     Sig(1:nx,j,k+1) = eta_final(k+1)*(Sig(1:nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                % end
                %
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                % end
                %
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                 
                Px(:,:,k+1) = zeros(nx);
                for jj = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))';
                end
                c = c + kappa*trace(Px(:,:,k+1));
                c_covariance = c_covariance + kappa*trace(Px(:,:,k+1));
                % for jj = 1:n_sampling_point
                % V_comp(:,:,k+1) = (Sig(1:nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))';
                % c = c + kappa*trace(w*V_comp(:,:,k+1));
                % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                
                % % debugging
                % c_variance(j,k+1) = kappa*trace(w*V_comp);
                %
                % V_comp_x = (Sig(1,j,k+1)-x_mean(1,k+1))*(Sig(1,j,k+1)-x_mean(1,k+1))';
                % c_variance_x(j,k+1) = kappa*trace(w*V_comp_x);
                % V_comp_xd = (Sig(7,j,k+1)-x_mean(7,k+1))*(Sig(7,j,k+1)-x_mean(7,k+1))';
                % c_variance_xd(j,k+1) = kappa*trace(w*V_comp_xd);
                %
                % V_comp_z = (Sig(3,j,k+1)-x_mean(3,k+1))*(Sig(3,j,k+1)-x_mean(3,k+1))';
                % c_variance_z(j,k+1) = kappa*trace(w*V_comp_z);
                % V_comp_zd = (Sig(9,j,k+1)-x_mean(9,k+1))*(Sig(9,j,k+1)-x_mean(9,k+1))';
                % c_variance_zd(j,k+1) = kappa*trace(w*V_comp_zd);
                % end
                
                % accumulate returned cost
                c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                c_quadratic = c_quadratic + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;
                %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                % % debugging
                % c_quadratic(k+1) = norm(x(:,k+1)-x_mean(:,k+1))^2;
                % c_quadratic_x(k+1) = norm(x(1,k+1)-x_mean(1,k+1))^2;
                % c_quadratic_xd(k+1) = norm(x(7,k+1)-x_mean(7,k+1))^2;
                % c_quadratic_z(k+1) = norm(x(3,k+1)-x_mean(3,k+1))^2;
                % c_quadratic_zd(k+1) = norm(x(7,k+1)-x_mean(9,k+1))^2;
                
                % derivative of variance matrix
                % gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
                dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                dmeanRdx(:,:,k+1) = zeros(obj.N,nx);
                dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                
                for j=k:-1:1
                    dTrVdx(j,:,k+1) = zeros(1,nx);
                    dTrVdu(j,:,k+1) = zeros(1,nu);
                    dmeanRdx(j,:,k+1) = zeros(1,nx);
                    dmeanRdu(j,:,k+1) = zeros(1,nu);
                    
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
                                    dSig_m_kplus1_dx = dfdx(:,:,m,2) + dfdSig(:,:,m,2);
                                end
                                dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]

                                while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                    chain_rule_indx = chain_rule_indx - 1;
                                end

                                % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
                                %     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
                                %     chain_rule_indx = chain_rule_indx - 1;
                                % end
                                %
                                % if (chain_rule_indx > 0)
                                %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                                %         chain_rule_indx = chain_rule_indx - 1;
                                %     end
                                % else
                                %     chain_rule_grad = eye(12);
                                % end
                                %
                                % dSig_m_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_dx;
                                % dSig_m_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_du;
                                
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
                            dSig_i_kplus1_dx = dfdx(:,:,i,2) + dfdSig(:,:,i,2);
                        end
                        dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                            dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                          dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                          dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
                          % % new sigma point due to resampling mechanism
                          % dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                          % dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                    end
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                    %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                end
            end
             
            dc = [];
            % cost gradient w.r.t x at first time step is zero
            for jj=1:obj.N % index for x_k
                dTrV_sum_dx_k = zeros(1, nx);
                if (jj == 1)
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(nx)-eye(nx));% equal to zero vector
                else
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                    %dmeanR_sum_dx_k = (inv(Px(:,:,jj))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                end
                                
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                    dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
                end
                dc = [dc, 1/2*dmeanR_sum_dx_k+dTrV_sum_dx_k];
                %dc = [dc, 1/2*dmeanR_sum_dx_k];
                %dc = [dc, dTrV_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            for jj=1:obj.N % index for u_k
                dTrV_sum_du_k = zeros(1, nu);
                dmeanR_sum_du_k = zeros(1, nu);
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                end
                dc = [dc, 1/2*dmeanR_sum_du_k+dTrV_sum_du_k];
                %dc = [dc, 1/2*dmeanR_sum_du_k];
                %dc = [dc, dTrV_sum_du_k];
            end
            
            % scale this robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
            fprintf('robust cost function: %4.8f\n',c);
            fprintf('robust mean deivation cost: %4.8f\n',obj.options.contact_robust_cost_coeff*c_quadratic);
            fprintf('robust covariance cost: %4.8f\n',obj.options.contact_robust_cost_coeff*c_covariance);
            tElapsed = toc(tStart);
            
            % figure(7),hold on;plot(c_quadratic_x,'b-')
            % figure(8),hold on;plot(c_quadratic_xd,'b-')
            % figure(9),hold on;plot(c_variance_x(1,:),'b-')
            % figure(10),hold on;plot(c_variance_xd(1,:),'b-')
            %
            % figure(11),hold on;plot(c_quadratic_z,'b-')
            % figure(12),hold on;plot(c_quadratic_zd,'b-')
            % figure(13),hold on;plot(c_variance_z(1,:),'b-')
            % figure(14),hold on;plot(c_variance_zd(1,:),'b-')
            
            % obj.cached_Px = Px;
            % fprintf('robust cost function: %4.8f\n',c);
            
            % % check gradient
            % disp('check gradient')
            % c_numeric = c;
            % dc_numeric = dc;
            %
            X0 = [x_full; Fext_full];
            % %X0 = X0 + randn(size(X0))*0.1;
            % 
            % fun = @(X0) robustVariancecost_check(obj, X0);
            % DerivCheck(fun, X0)
            % disp('finish')
            c_numeric = c;
            dc_numeric = dc;
             
            % [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_scaled_check(obj,X0),X0,struct('grad_method','numerical'));
            %
            % [c_numeric,dc_numeric] = robustVariancecost_check(obj, X0);
            %
            % valuecheck(dc,dc_numeric,1e-5);
            % valuecheck(c,c_numeric,1e-5);
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %
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
                dd =  dot(rr, JJ)-f2+f1;
                fprintf('difference between numerical and analytical: %4.15f\n',dd);
            end
            
            function [c,dc] = robustVariancecost_scaled_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                
                manip = obj.plant.getManipulator();
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nx = nq+nv;
                nu = obj.nFext;%obj.plant.getNumInputs;
                obj.nu = nu;
                
                % sigma points
                Px = zeros(nx,nx,obj.N);
                Px(:,:,1) = obj.cached_Px(:,:,1);
                Px_init = Px(:,:,1);
                
                if strcmp(manip.uncertainty_source,'friction_coeff')
                    w_mu = manip.uncertain_mu_set;
                    w_noise = [w_mu];
                    Pw = diag([0.01]);
                elseif strcmp(manip.uncertainty_source,'terrain_height')
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_phi];
                    Pw = diag([0.016]);
                elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                    w_mu = manip.uncertain_mu_set;
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_mu;w_phi];
                    Pw = diag([0.01, 0.016]);
                elseif isempty(manip.uncertainty_source)
                    Pw = [];
                    w_noise = [];
                elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                    %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                    w_noise = [w_mu];
                end
                
                % disturbance variance
                % currently only consider terrain height and/or friction coefficient
                nw = size(Pw,1);
                sampling_method = 2;%option 1: unscented transform, option 2: random sampling with a smaller number
                if sampling_method == 1
                    n_sampling_point = 2*(nx+nw);
                    scale = .01;% [to be tuned]
                    w = 0.5/scale^2;
                elseif sampling_method == 2
                    n_sampling_point = 5;
                    w_state = load('state_noise.dat');%std:0.001
                    %w_state = 0.001*randn(28,62);
                    %save -ascii state_noise_small.dat w_state
                    w = 1;
                end
                w_avg = 1/n_sampling_point;
                K = obj.options.K;
                alpha = obj.options.alpha;
                
                %initialize c and dc
                kappa = obj.options.kappa;
                x_mean = zeros(nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
                c = 0;
                %c = kappa*trace(Px_init);
                %c_quadratic = 0;
                %c_variance = 0;
                dc = zeros(1, 1+obj.N*(nx+nu));
                
                % initialize gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,1) = zeros(obj.N-1,nx);
                dTrVdu(:,:,1) = zeros(obj.N-1,nu);
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(timestep_updated,Sig,u_fdb_k)
                    [xdn,df] = obj.plant.update(timestep_updated,Sig,u_fdb_k);
                end
                
                plant_update = @objPlantUpdate;
                obj.plant.time_step = timestep_updated;
                df_previous_full = [];
                xdn_previous_full = [];
                noise_sample_type = 1;
                
                for k = 1:obj.N-1%[Ye: double check the index]
                    if sampling_method == 1
                        %Propagate sigma points through nonlinear dynamics
                        if k == 1
                            [S,d] = chol(blkdiag(Px_init, Pw), 'lower');
                            
                            if d
                                diverge = k;
                                return;
                            end
                            S = scale*S;
                            Sig_init(:,:,k) = [S -S];
                            
                            %state dimensionality of sigma point depends on the number of uncertainties.
                            %Our case w_noise has 0, 1, or 2 dimensions
                            if isempty(w_noise)
                                for j = 1:n_sampling_point
                                    Sig_init(:,j,k) = Sig_init(:,j,k) + x(:,k);
                                end
                            else
                                for j = 1:n_sampling_point
                                    Sig_init(:,j,k) = Sig_init(:,j,k) + [x(:,k); w_noise(:,k)];
                                end
                            end
                            Sig(:,:,k) = Sig_init(:,:,k);
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
                            
                            Px_init = zeros(nx);
                            for j = 1:n_sampling_point
                                Px_init = Px_init + w*(Sig(1:nx,j,k)-x_mean(:,k))*(Sig(1:nx,j,k)-x_mean(:,k))';
                            end
                            c = c + trace(Px_init);
                        else
                            for j = 1:n_sampling_point
                                Sig_init(:,j,k) = (1-alpha)*Sig(:,j,k) + alpha*x(:,k);
                            end
                        end
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:n_sampling_point
                        if strcmp(manip.uncertainty_source, 'friction_coeff')
                            uncertain_mu = w_mu(j);
                        elseif strcmp(manip.uncertainty_source, 'terrain_height')
                            terrain_index = j;
                        elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                            uncertain_mu = w_mu(j);
                            terrain_index = j;
                        end
                        
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig_init(1:nx/2,j,k),Sig_init(nx/2+1:nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        
                        % add feedback control
                        t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                        u_fdb_k = u(:,k) - K*(Sig_init(1:nx,j,k) - x(:,k));
                        
                        if k > 1
                            x_previous = xdn_previous_full(:,j);
                            df_previous = df_previous_full(:,:,j);
                        end
                        
                        if noise_sample_type == 1
                            [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                        elseif noise_sample_type == 2
                            xdn_analytical(:,j) = zeros(nx,1);
                            df_analytical(:,:,j) = zeros(nx,1+nx+nu);
                            for kk=1:length(w_noise)
                                [xdn_analytical_sample(:,j),df_analytical_sample(:,:,j)] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k);
                                xdn_analytical(:,j) = xdn_analytical(:,j) + xdn_analytical_sample(:,j);
                                df_analytical(:,:,j) = df_analytical(:,:,j) + df_analytical_sample(:,:,j);
                            end
                            xdn_analytical(:,j) = xdn_analytical(:,j)/length(w_noise);
                            df_analytical(:,:,j) = df_analytical(:,:,j)/length(w_noise);
                        end
                        % %numerical diff
                        % dt = diag(max(sqrt(eps(timestep_updated)), 1e-7));
                        % dx = diag(max(sqrt(eps(Sig_init(1:nx,j,k))), 1e-7*ones(nx,1)));
                        % du = diag(max(sqrt(eps(u_fdb_k)), 1e-7*ones(nu,1)));
                        %
                        % [xdnp,~] = feval(plant_update,timestep_updated+dt,Sig_init(1:nx,j,k),u_fdb_k);
                        % [xdnm,~] = feval(plant_update,timestep_updated-dt,Sig_init(1:nx,j,k),u_fdb_k);
                        % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                        %
                        % N_finite_diff_x = length(Sig_init(1:nx,j,k));
                        % for m = 1:N_finite_diff_x
                        %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k)+dx(:,m),u_fdb_k);
                        %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k)-dx(:,m),u_fdb_k);
                        %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                        % end
                        %
                        % N_finite_diff_u = length(u_fdb_k);
                        % for m = 1:N_finite_diff_u
                        %     [xdnp,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k+du(:,m));
                        %     [xdnm,~] = feval(plant_update,timestep_updated,Sig_init(1:nx,j,k),u_fdb_k-du(:,m));
                        %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                        % end
                        
                        % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                        %     keyboard
                        % end
                        
                        xdn(:,j) = xdn_analytical(:,j);
                        df(:,:,j) = df_analytical(:,:,j);
                        
                        Sig(1:nx,j,k+1) = xdn(1:nx,j);
                        dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                        dfdSig(:,:,j,k+1) = (1-alpha)*df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*(1-alpha)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*(1-alpha)*K + alpha*df(:,2:nx+1,j);
                    end
                    xdn_previous_full = xdn;
                    df_previous_full = df;
                    
                    %calculate mean and variance w.r.t. [x_k] from sigma points
                    x_mean(:,k+1) = zeros(nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                    end
                    
                    % % shift sigma point
                    % for j = 1:n_sampling_point
                    %     Sig(1:nx,j,k+1) = Sig(1:nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                    %     eta_optimal = 1;%etaEstimation(Sig(1:nx,j,k+1),x(:,k+1));
                    %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                    %         eta_final(k+1) = eta_optimal;
                    %     end
                    % end
                    %
                    % % scale Sigma point deviation by eta
                    % for j = 1:n_sampling_point
                    %     Sig(1:nx,j,k+1) = eta_final(k+1)*(Sig(1:nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                    % end
                    %
                    % % recalculate mean and variance w.r.t. [x_k] from sigma points
                    % x_mean(:,k+1) = zeros(nx,1);
                    % for j = 1:n_sampling_point
                    %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:nx,j,k+1);
                    % end
                    %
                    % % check that the mean deviation term is cancelled out
                    % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                    %     disp('shifting scheme is not correct')
                    %     keyboard
                    % end
                    
                    Px(:,:,k+1) = zeros(nx);
                    for jj = 1:n_sampling_point
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))';
                    end
                    c = c + kappa*trace(Px(:,:,k+1));
                    
                    % for jj = 1:n_sampling_point
                    % V_comp(:,:,k+1) = (Sig(1:nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:nx,jj,k+1)-x_mean(:,k+1))';
                    % c = c + kappa*trace(w*V_comp(:,:,k+1));
                    % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                    
                    % % debugging
                    % c_variance(j,k+1) = kappa*trace(w*V_comp);
                    %
                    % V_comp_x = (Sig(1,j,k+1)-x_mean(1,k+1))*(Sig(1,j,k+1)-x_mean(1,k+1))';
                    % c_variance_x(j,k+1) = kappa*trace(w*V_comp_x);
                    % V_comp_xd = (Sig(7,j,k+1)-x_mean(7,k+1))*(Sig(7,j,k+1)-x_mean(7,k+1))';
                    % c_variance_xd(j,k+1) = kappa*trace(w*V_comp_xd);
                    %
                    % V_comp_z = (Sig(3,j,k+1)-x_mean(3,k+1))*(Sig(3,j,k+1)-x_mean(3,k+1))';
                    % c_variance_z(j,k+1) = kappa*trace(w*V_comp_z);
                    % V_comp_zd = (Sig(9,j,k+1)-x_mean(9,k+1))*(Sig(9,j,k+1)-x_mean(9,k+1))';
                    % c_variance_zd(j,k+1) = kappa*trace(w*V_comp_zd);
                    % end
                    
                    % accumulate returned cost
                    c = c + 1/2*norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                    % % debugging
                    % c_quadratic(k+1) = norm(x(:,k+1)-x_mean(:,k+1))^2;
                    % c_quadratic_x(k+1) = norm(x(1,k+1)-x_mean(1,k+1))^2;
                    % c_quadratic_xd(k+1) = norm(x(7,k+1)-x_mean(7,k+1))^2;
                    % c_quadratic_z(k+1) = norm(x(3,k+1)-x_mean(3,k+1))^2;
                    % c_quadratic_zd(k+1) = norm(x(7,k+1)-x_mean(9,k+1))^2;
                    
                    % derivative of variance matrix
                    % gradient of Tr(V) w.r.t state vector x
                    dTrVdx(:,:,k+1) = zeros(obj.N-1,nx);
                    dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
                    dmeanRdx(:,:,k+1) = zeros(obj.N,nx);
                    dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
                    
                    for j=k:-1:1
                        dTrVdx(j,:,k+1) = zeros(1,nx);
                        dTrVdu(j,:,k+1) = zeros(1,nu);
                        dmeanRdx(j,:,k+1) = zeros(1,nx);
                        dmeanRdu(j,:,k+1) = zeros(1,nu);
                        
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
                                        dSig_m_kplus1_dx = dfdx(:,:,m,2) + dfdSig(:,:,m,2);
                                    end
                                    dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                    
                                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                        dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                        dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                        chain_rule_indx = chain_rule_indx - 1;
                                    end
                                    
                                    % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    %     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
                                    %     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
                                    %     chain_rule_indx = chain_rule_indx - 1;
                                    % end
                                    %
                                    % if (chain_rule_indx > 0)
                                    %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                                    %         chain_rule_indx = chain_rule_indx - 1;
                                    %     end
                                    % else
                                    %     chain_rule_grad = eye(12);
                                    % end
                                    %
                                    % dSig_m_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_dx;
                                    % dSig_m_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_du;
                                    
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
                                dSig_i_kplus1_dx = dfdx(:,:,i,2) + dfdSig(:,:,i,2);
                            end
                            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            
                            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
                            % % new sigma point due to resampling mechanism
                            % dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            % dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        end
                        
                        % gradient of mean residual w.r.t state x and control u, assume norm 2
                        dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                        dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                        %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                        %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    end
                end
                
                dc = [];
                % cost gradient w.r.t x at first time step is zero
                for jj=1:obj.N % index for x_k
                    dTrV_sum_dx_k = zeros(1, nx);
                    if (jj == 1)
                        dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(nx)-eye(nx));% equal to zero vector
                    else
                        dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                        %dmeanR_sum_dx_k = (inv(Px(:,:,jj))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                    end
                    
                    for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                        dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
                    end
                    dc = [dc, 1/2*dmeanR_sum_dx_k+dTrV_sum_dx_k];
                    %dc = [dc, 1/2*dmeanR_sum_dx_k];
                    %dc = [dc, dTrV_sum_dx_k];
                end
                
                % cost gradient w.r.t u at first time step is zero, since
                for jj=1:obj.N % index for u_k
                    dTrV_sum_du_k = zeros(1, nu);
                    dmeanR_sum_du_k = zeros(1, nu);
                    for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                        dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                    end
                    dc = [dc, 1/2*dmeanR_sum_du_k+dTrV_sum_du_k];
                    %dc = [dc, 1/2*dmeanR_sum_du_k];
                    %dc = [dc, dTrV_sum_du_k];
                end
                
                % scale this robust cost
                c = obj.options.contact_robust_cost_coeff*c;
                dc = obj.options.contact_robust_cost_coeff*dc;
            end
        end
        
        function [c,dc] = robustVariancecost_ML_multi_step_prop(obj, x_full, Fext_full)
            global timestep_updated
            global time_step
            global x_previous
            global df_previous
            global terrain_index
            global uncertain_mu
            manip = obj.plant.getManipulator();
            
            x = reshape(x_full, obj.nx, obj.N);
            u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
            nq = obj.plant.getNumPositions;
            nv = obj.plant.getNumVelocities;
            nx = nq+nv;
            nu = obj.nFext;%obj.plant.getNumInputs;
            obj.nu = nu;
            
            % sigma points
            Px = zeros(obj.nx,obj.nx,obj.N);
            Px(:,:,1) = obj.cached_Px(:,:,1);
            
            if strcmp(manip.uncertainty_source,'friction_coeff')
                w_mu = manip.uncertain_mu_set;
                w_noise = [w_mu];
                Pw = diag([0.01]);
            elseif strcmp(manip.uncertainty_source,'terrain_height')
                w_phi = manip.uncertain_position_set;
                w_noise = [w_phi];
                Pw = diag([0.016]);
            elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                w_mu = manip.uncertain_mu_set;
                w_phi = manip.uncertain_position_set;
                w_noise = [w_mu;w_phi];
                Pw = diag([0.01, 0.016]);
            elseif isempty(manip.uncertainty_source)
                Pw = [];
                w_noise = [];
            elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                save -ascii friction_coeff_noise.dat w_mu
                %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                %w_phi = [x;y];%height noise
                %save -ascii initial_position_noise.dat w_phi
                w_noise = [w_mu];
            end
            
            % disturbance variance
            % currently only consider terrain height and/or friction coefficient
            scale = .01;% [to be tuned]
            w = 0.5/scale^2;
            nw = size(Pw,1);
            n_sampling_point = 2*(obj.nx+nw);
            w_avg = 1/n_sampling_point;
            
            K = obj.options.K;
            eta_final = ones(obj.N,1);
            
            %initialize c and dc
            kappa = obj.options.kappa;
            x_mean = zeros(obj.nx, obj.N);
            % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
            c = 0;
            %c = kappa*trace(Px(:,:,1));
            %c_quadratic = 0;
            %c_variance = 0;
            dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
            
            % initialize gradient of Tr(V) w.r.t state vector x
            dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
            dTrVdu(:,:,1) = zeros(obj.N-1,nu);
            
            % time counter
            tStart = tic;
            
            function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
            end
            
            plant_update = @objPlantUpdate;
            %obj.plant.time_step = time_step;
            df_previous_full = [];
            xdn_previous_full = [];
            
            for k = 1:obj.N-1%[Ye: double check the index]
                %Propagate sigma points through nonlinear dynamics
                if k == 1
                    [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                    if d
                        diverge = k;
                        return;
                    end
                    
                    S = scale*S;
                    Sig(:,:,k) = [S -S];
                    if isempty(w_noise)
                        for j = 1:n_sampling_point
                            Sig(:,j,k) = Sig(:,j,k) + [x(:,k)];
                        end
                    else
                        for j = 1:n_sampling_point
                            Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                        end
                    end
                    x_mean(:,k) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:obj.nx,j,k);
                    end
                    %c = c + norm(x(:,k)-x_mean(:,k))^2;
                    c = c + 1/2*log(det(Px(:,:,k)));
                end
                
                %Propagate sigma points through nonlinear dynamics
                for j = 1:n_sampling_point
                    noise_index = j;
                    
                    % a hacky way to implement the control input
                    [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                    Hinv(:,:,j,k) = inv(H);
                    
                    if strcmp(manip.uncertainty_source, 'friction_coeff')
                        uncertain_mu = w_mu(j);
                    elseif strcmp(manip.uncertainty_source, 'terrain_height')
                        terrain_index = j;
                    elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                        uncertain_mu = w_mu(j);
                        terrain_index = j;
                    end
                    
                    % add feedback control
                    t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                    u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                    
                    if k > 1
                        x_previous = xdn_previous_full(:,j);
                        df_previous = df_previous_full(:,:,j);
                    end
                    
                    % estimate whether current state is close to contact
                    [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                    phi_bottom = phi_current(2:2:end);
                    active_threshold = 0.5;
                    if 0 %any(phi_bottom<active_threshold)
                        for kk = 1:size_terrain_sample
                            if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                obj.plant.friction_coeff = w_mu(kk);
                            elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                obj.plant.terrain_index = kk;
                            end
                            [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                        end
                        xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                        Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                        Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                        
                    else
                        [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig(1:nx,j,k),u_fdb_k);
                        
                        % %numerical diff
                        % dt = diag(sqrt(eps(t)));
                        % dx = diag(max(sqrt(eps(Sig(1:obj.nx,j,k))), 1e-3*ones(obj.nx,1)));
                        % du = diag(max(sqrt(eps(u_fdb_k)), 1e-3*ones(nu,1)));
                        %
                        % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                        % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                        %
                        % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                        % for m = 1:N_finite_diff_x
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                        %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                        % end
                        %
                        % N_finite_diff_u = length(u_fdb_k);
                        % for m = 1:N_finite_diff_u
                        %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                        %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                        %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                        % end
                        %
                        % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 1e-2)
                        %     keyboard
                        % end
                        
                        xdn = xdn_analytical;
                        df(:,:,j) = df_analytical(:,:,j);
                        
                        Sig(1:nx,j,k+1) = xdn(1:nx,j);
                        dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                        dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                    end
                end
                xdn_previous_full = xdn;
                df_previous_full = df;
                
                %calculate mean and variance w.r.t. [x_k] from sigma points
                x_mean(:,k+1) = zeros(obj.nx,1);
                for j = 1:n_sampling_point
                    x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                end
                
                % % shift sigma point
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                %         eta_final(k+1) = eta_optimal;
                %     end
                % end
                %
                % % scale Sigma point deviation by eta
                % for j = 1:n_sampling_point
                %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                % end
                %
                % % recalculate mean and variance w.r.t. [x_k] from sigma points
                % x_mean(:,k+1) = zeros(obj.nx,1);
                % for j = 1:n_sampling_point
                %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                % end
                %
                % % check that the mean deviation term is cancelled out
                % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                %     disp('shifting scheme is not correct')
                %     keyboard
                % end
                
                Px(:,:,k+1) = zeros(obj.nx);
                for jj = 1:n_sampling_point
                    Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                end
                %c = c + kappa*trace(Px(:,:,k+1));
                
                % accumulate returned cost
                %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                
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
                    
                    tic
                    % gradient w.r.t state x
                    dSig_m_kplus1_dx_sum = zeros(obj.nx);
                    % gradient w.r.t control u
                    dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                    dSig_i_kplus1_dx_resample = zeros(obj.nx,obj.nx,n_sampling_point);
                    dSig_i_kplus1_du_resample = zeros(obj.nx,nu,n_sampling_point);
                    
                    for i=1:n_sampling_point
                        if i == 1
                            for m = 1:n_sampling_point% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                % gradient of Tr(V_{k+1}) w.r.t control x and u
                                dSig_m_kplus1_dx = zeros(obj.nx);
                                dSig_m_kplus1_du = zeros(obj.nx,1);
                                
                                chain_rule_indx = k-j;
                                if j ~= 1
                                    %dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                    dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                                else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                    dSig_m_kplus1_dx = dfdx(:,:,m,2) + dfdSig(:,:,m,2);
                                end
                                dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]

                                while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                    dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                                    chain_rule_indx = chain_rule_indx - 1;
                                end

                                % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %     dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
                                %     dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
                                %     chain_rule_indx = chain_rule_indx - 1;
                                % end
                                %
                                % if (chain_rule_indx > 0)
                                %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                                %         chain_rule_indx = chain_rule_indx - 1;
                                %     end
                                % else
                                %     chain_rule_grad = eye(12);
                                % end
                                %
                                % dSig_m_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_dx;
                                % dSig_m_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_m_kplus1_du;
                                
                                dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                            end
                        end
                        
                        % run 2*(obj.nx+nw) times in total to obtain
                        % gradient w.r.t sigma points
                        dSig_i_kplus1_dx = zeros(obj.nx);
                        dSig_i_kplus1_du = zeros(obj.nx,nu);
                        chain_rule_indx = k-j;
                        if j ~= 1
                            %dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,j+1);
                            dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
                        else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                            dSig_i_kplus1_dx = dfdx(:,:,i,2) + dfdSig(:,:,i,2)*eye(obj.nx);
                        end
                        dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                        
                        while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                            dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                            dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                            chain_rule_indx = chain_rule_indx - 1;
                        end
                        
                        % while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        %     dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_dx;
                        %     dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_i_kplus1_du;
                        %     chain_rule_indx = chain_rule_indx - 1;
                        % end
                        %
                        % if (chain_rule_indx > 0)
                        %     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        %         chain_rule_grad = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*chain_rule_grad;
                        %         chain_rule_indx = chain_rule_indx - 1;
                        %     end
                        % else
                        %     chain_rule_grad = eye(obj.nx);
                        % end
                        %
                        % dSig_i_kplus1_dx = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_dx;
                        % dSig_i_kplus1_du = chain_rule_grad*(1-w_avg)*dSig_i_kplus1_du;
                        
                        dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                        dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        
                        % % new sigma point due to resampling mechanism
                        % dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                        % dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                    end
                    
                    % %% resample scheme
                    % dSig_kplus1_dx_sum_resample = zeros(obj.nx,obj.nx);
                    % dSig_kplus1_du_sum_resample = zeros(obj.nx,nu);
                    %
                    % for i =1:n_sampling_point
                    %     dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
                    %     dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
                    % end
                    %
                    % for i =1:n_sampling_point
                    %     dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
                    %     dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
                    % end
                    % %% resample scheme
                    
                    % gradient of mean residual w.r.t state x and control u, assume norm 2
                    dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                    dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                    %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                    %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                    
                    % dCovdx = zeros(nx,nx,nx);
                    % dCovdu = zeros(nx,nx,nu);
                    % toc
                    %
                    % tic
                    % for pp=1:nx
                    %     for i=1:n_sampling_point
                    %         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                    %         for nn=1:nx
                    %             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                    %             if nn == 1
                    %                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                    %                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                    %             elseif nn == nx
                    %                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                    %             else
                    %                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                    %                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                    %                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                    %             end
                    %             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_resample(nn,pp,i) - w_avg*dSig_kplus1_dx_sum_resample(nn,pp));
                    %             if pp <= nu
                    %                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_resample(nn,pp,i) - w_avg*dSig_kplus1_du_sum_resample(nn,pp));
                    %             end
                    %         end
                    %     end
                    %     dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                    %     if pp <= nu
                    %         dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                    %     end
                    % end
                    % toc
                end
            end
            
            dc = [];
            % cost gradient w.r.t x at first time step is zero
            for jj=1:obj.N % index for x_k
                dTrV_sum_dx_k = zeros(1, obj.nx);
                if (jj == 1)
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                else
                    dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))';%i.i.d mean deviation version
                    %dmeanR_sum_dx_k = (pinv(Px(:,:,jj))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                end
                
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                    dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
                end
                %dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                dc = [dc, dmeanR_sum_dx_k];
                %dc = [dc, kappa*dTrV_sum_dx_k];
            end
            
            % cost gradient w.r.t u at first time step is zero, since
            for jj=1:obj.N % index for u_k
                dTrV_sum_du_k = zeros(1, nu);
                dmeanR_sum_du_k = zeros(1, nu);
                for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                    dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                    dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                end
                %dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                dc = [dc, dmeanR_sum_du_k];
                %dc = [dc, kappa*dTrV_sum_du_k];
            end
            
            % scale this robust cost
            c = obj.options.contact_robust_cost_coeff*c;
            dc = obj.options.contact_robust_cost_coeff*dc;
            fprintf('robust cost function: %4.8f\n',c);
            tElapsed = toc(tStart);
            
            % figure(7),hold on;plot(c_quadratic_x,'b-')
            % figure(8),hold on;plot(c_quadratic_xd,'b-')
            % figure(9),hold on;plot(c_variance_x(1,:),'b-')
            % figure(10),hold on;plot(c_variance_xd(1,:),'b-')
            %
            % figure(11),hold on;plot(c_quadratic_z,'b-')
            % figure(12),hold on;plot(c_quadratic_zd,'b-')
            % figure(13),hold on;plot(c_variance_z(1,:),'b-')
            % figure(14),hold on;plot(c_variance_zd(1,:),'b-')
            %
            % figure(15)
            % clf
            % Sig_permute = permute(Sig,[1,3,2]);
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(1,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(1,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,7])
            %
            % figure(16)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(7,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(7,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,11])
            %
            % figure(17)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(3,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(3,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([0,4])
            %
            % figure(18)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(9,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(9,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([-7,0])
            %
            % figure(19)
            % clf
            % for j = 1:n_sampling_point
            %     plot(Sig_permute(5,:,j),'r-')
            %     hold on;
            % end
            % hold on;
            % plot(x(5,:),'b-','Linewidth',3)
            % xlim([0,30]);ylim([-1,1])
            %
            % obj.cached_Px = Px;
            % fprintf('robust cost function: %4.8f\n',c);
            
            % % check gradient
            % disp('check gradient')
            % c_numeric = c;
            % dc_numeric = dc;
            %
            X0 = [x_full; Fext_full];
            % %X0 = X0 + randn(size(X0))*0.1;
            %
            % fun = @(X0) robustVariancecost_check(obj, X0);
            % DerivCheck(fun, X0)
            % disp('finish')
            c_numeric = c;
            dc_numeric = dc;
            
            % [c_numeric,dc_numeric] = geval(@(X0) robustVariancecost_check(obj,X0),X0,struct('grad_method','numerical'));
            %
            % [c_numeric,dc_numeric] = robustVariancecost_check(obj, X0);
            %
            % valuecheck(dc,dc_numeric,1e-5);
            % valuecheck(c,c_numeric,1e-5);
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %
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
                dd =  dot(rr, JJ)-f2+f1;
                fprintf('difference between numerical and analytical: %4.15f\n',dd);
            end
            
            function [c,dc] = robustVariancecost_check(obj, X0)
                x_full = X0(1:obj.nx*obj.N);
                Fext_full = X0(obj.nx*obj.N+1:end);
                
                manip = obj.plant.getManipulator();
                
                x = reshape(x_full, obj.nx, obj.N);
                u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
                nq = obj.plant.getNumPositions;
                nv = obj.plant.getNumVelocities;
                nx = nq+nv;
                nu = obj.nFext;%obj.plant.getNumInputs;
                obj.nu = nu;
                
                % sigma points
                Px = zeros(obj.nx,obj.nx,obj.N);
                Px(:,:,1) = obj.cached_Px(:,:,1);
                
                if strcmp(manip.uncertainty_source,'friction_coeff')
                    w_mu = manip.uncertain_mu_set;
                    w_noise = [w_mu];
                    Pw = diag([0.01]);
                elseif strcmp(manip.uncertainty_source,'terrain_height')
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_phi];
                    Pw = diag([0.016]);
                elseif strcmp(manip.uncertainty_source,'friction_coeff+terrain_height')
                    w_mu = manip.uncertain_mu_set;
                    w_phi = manip.uncertain_position_set;
                    w_noise = [w_mu;w_phi];
                    Pw = diag([0.01, 0.016]);
                elseif isempty(manip.uncertainty_source)
                    Pw = [];
                    w_noise = [];
                elseif strcmp(manip.uncertainty_source,'generate_new_noise_set')
                    w_mu = normrnd(ones(1,n_sampling_point),sqrt(Pw(1,1)),1,n_sampling_point);
                    save -ascii friction_coeff_noise.dat w_mu
                    %x = (1-2*rand(1,n_sampling_point))*sqrt(Pw(1,1));
                    %y = (1-2*rand(1,n_sampling_point))*sqrt(Pw(2,2));
                    %w_phi = [x;y];%height noise
                    %save -ascii initial_position_noise.dat w_phi
                    w_noise = [w_mu];
                end
                
                % disturbance variance
                % currently only consider terrain height and/or friction coefficient
                scale = .01;% [to be tuned]
                w = 0.5/scale^2;
                nw = size(Pw,1);
                n_sampling_point = 2*(obj.nx+nw);
                w_avg = 1/n_sampling_point;
                
                K = obj.options.K;
                eta_final = ones(obj.N,1);
                
                %initialize c and dc
                kappa = obj.options.kappa;
                x_mean = zeros(obj.nx, obj.N);
                % mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
                %c = 0;
                c = kappa*trace(Px(:,:,1));
                %c_quadratic = 0;
                %c_variance = 0;
                dc = zeros(1, 1+obj.N*(obj.nx+obj.nu));
                
                % initialize gradient of Tr(V) w.r.t state vector x
                dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
                dTrVdu(:,:,1) = zeros(obj.N-1,nu);
                
                % time counter
                tStart = tic;
                
                function [xdn,df] = objPlantUpdate(noise_index,Sig,u_fdb_k)
                    [xdn,df] = obj.plant.update(noise_index,Sig,u_fdb_k);
                end
                
                plant_update = @objPlantUpdate;
                %obj.plant.time_step = time_step;
                df_previous_full = [];
                xdn_previous_full = [];
                
                for k = 1:obj.N-1%[Ye: double check the index]
                    %Propagate sigma points through nonlinear dynamics
                    if k == 1
                        [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
                        if d
                            diverge = k;
                            return;
                        end
                        S = scale*S;
                        Sig(:,:,k) = [S -S];
                        if isempty(w_noise)
                            for j = 1:n_sampling_point
                                Sig(:,j,k) = Sig(:,j,k) + [x(:,k)];
                            end
                        else
                            for j = 1:n_sampling_point
                                Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
                            end
                        end
                        x_mean(:,k) = zeros(obj.nx,1);
                        for j = 1:n_sampling_point
                            x_mean(:,k) = x_mean(:,k) + w_avg*Sig(1:obj.nx,j,k);
                        end
                        %c = c + norm(x(:,k)-x_mean(:,k))^2;
                        
                        % % debugging
                        % c_quadratic(k) = norm(x(:,k)-x_mean(:,k))^2;
                        % c_quadratic_x(k) = norm(x(1,k)-x_mean(1,k))^2;
                        % c_quadratic_xd(k) = norm(x(7,k)-x_mean(7,k))^2;
                        % c_quadratic_z(k) = norm(x(3,k)-x_mean(3,k))^2;
                        % c_quadratic_zd(k) = norm(x(9,k)-x_mean(9,k))^2;
                        % for j = 1:n_sampling_point
                        %     V_comp = (Sig(1:obj.nx,j,k)-x_mean(:,k))*(Sig(1:obj.nx,j,k)-x_mean(:,k))';
                        %     c_variance(j,k) = kappa*trace(w*V_comp);
                        %
                        %     % debugging
                        %     V_comp_x = (Sig(1,j,k)-x_mean(1,k))*(Sig(1,j,k)-x_mean(1,k))';
                        %     c_variance_x(j,k) = kappa*trace(w*V_comp_x);
                        %     V_comp_xd = (Sig(7,j,k)-x_mean(7,k))*(Sig(7,j,k)-x_mean(7,k))';
                        %     c_variance_xd(j,k) = kappa*trace(w*V_comp_xd);
                        %
                        %     V_comp_z = (Sig(3,j,k)-x_mean(3,k))*(Sig(3,j,k)-x_mean(3,k))';
                        %     c_variance_z(j,k) = kappa*trace(w*V_comp_z);
                        %     V_comp_zd = (Sig(9,j,k)-x_mean(9,k))*(Sig(9,j,k)-x_mean(9,k))';
                        %     c_variance_zd(j,k) = kappa*trace(w*V_comp_zd);
                        % end
                    end
                    
                    %Propagate sigma points through nonlinear dynamics
                    for j = 1:n_sampling_point
                        noise_index = j;
                        
                        % a hacky way to implement the control input
                        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
                        Hinv(:,:,j,k) = inv(H);
                        
                        if strcmp(manip.uncertainty_source, 'friction_coeff')
                            uncertain_mu = w_mu(j);
                        elseif strcmp(manip.uncertainty_source, 'terrain_height')
                            terrain_index = j;
                        elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                            uncertain_mu = w_mu(j);
                            terrain_index = j;
                        end
                        
                        % add feedback control
                        t = timestep_updated*(k-1);%[double make sure obj.h is updated correctly]
                        u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
                        
                        if k > 1
                            x_previous = xdn_previous_full(:,j);
                            df_previous = df_previous_full(:,:,j);
                        end
                        
                        % estimate whether current state is close to contact
                        [phi_current,~,~,~,~,~,~,~,~,~,~,~] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j,k),false,obj.options.active_collision_options);
                        phi_bottom = phi_current(2:2:end);
                        active_threshold = 0.5;
                        if 0 %any(phi_bottom<active_threshold)
                            for kk = 1:size_terrain_sample
                                if strcmp(obj.plant.uncertainty_source, 'friction_coeff')
                                    obj.plant.friction_coeff = w_mu(kk);
                                elseif strcmp(obj.plant.uncertainty_source, 'terrain_height')
                                    obj.plant.terrain_index = kk;
                                end
                                [xdn_sample(:,kk),df_sample(:,:,kk)] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
                            end
                            xdn_mean = sum(xdn_sample,2)/size_terrain_sample;
                            Sig(1:obj.nx/2,j,k+1) = xdn_mean(1:obj.nx/2);
                            Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn_mean(obj.nx/2+1:obj.nx);
                            
                        else
                            [xdn_analytical(:,j),df_analytical(:,:,j)] = feval(plant_update,noise_index,Sig(1:nx,j,k),u_fdb_k);
                            
                            % %numerical diff
                            % dt = diag(sqrt(eps(t)));
                            % dx = diag(sqrt(eps(Sig(1:obj.nx,j,k))));
                            % du = diag(sqrt(eps(u_fdb_k)));
                            %
                            % [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k);
                            % df_numeric(:,1) = (xdnp-xdnm)/(2*dt);
                            %
                            % N_finite_diff_x = length(Sig(1:obj.nx,j,k));
                            % tic
                            % for m = 1:N_finite_diff_x
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)+dx(:,m),u_fdb_k);
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k)-dx(:,m),u_fdb_k);
                            %     df_numeric(:,m+1) = (xdnp-xdnm)/(2*dx(m,m));
                            % end
                            %
                            % N_finite_diff_u = length(u_fdb_k);
                            %
                            % for m = 1:N_finite_diff_u
                            %     [xdnp,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k+du(:,m));
                            %     [xdnm,~] = feval(plant_update,noise_index,Sig(1:obj.nx,j,k),u_fdb_k-du(:,m));
                            %     df_numeric(:,m+1+N_finite_diff_x) = (xdnp-xdnm)/(2*du(m,m));
                            % end
                            %
                            % if (sum(sum(abs(df_analytical(:,:,j) - df_numeric))) > 10)
                            %     keyboard
                            % end
                            %
                            xdn = xdn_analytical;
                            df(:,:,j) = df_analytical(:,:,j);
                            
                            Sig(1:nx,j,k+1) = xdn(1:nx,j);
                            dfdu(:,:,j,k+1) = df(:,end-nu+1:end,j);
                            dfdSig(:,:,j,k+1) = df(:,2:nx+1,j) - dfdu(:,:,j,k+1)*K;
                            dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
                        end
                    end
                    xdn_previous_full = xdn;
                    df_previous_full = df;
                    
                    %calculate mean and variance w.r.t. [x_k] from sigma points
                    x_mean(:,k+1) = zeros(obj.nx,1);
                    for j = 1:n_sampling_point
                        x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    end
                    
                    % % shift sigma point
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = Sig(1:obj.nx,j,k+1) - x_mean(:,k+1) + x(:,k+1);
                    %     eta_optimal = 1;%etaEstimation(Sig(1:obj.nx,j,k+1),x(:,k+1));
                    %     if eta_optimal < eta_final(k+1)% actually represent for k^th time step, use k+1 for notation consistency later
                    %         eta_final(k+1) = eta_optimal;
                    %     end
                    % end
                    %
                    % % scale Sigma point deviation by eta
                    % for j = 1:n_sampling_point
                    %     Sig(1:obj.nx,j,k+1) = eta_final(k+1)*(Sig(1:obj.nx,j,k+1) - x(:,k+1)) + x(:,k+1);
                    % end
                    %
                    % % recalculate mean and variance w.r.t. [x_k] from sigma points
                    % x_mean(:,k+1) = zeros(obj.nx,1);
                    % for j = 1:n_sampling_point
                    %     x_mean(:,k+1) = x_mean(:,k+1) + w_avg*Sig(1:obj.nx,j,k+1);
                    % end
                    %
                    % % check that the mean deviation term is cancelled out
                    % if any(abs(x_mean(:,k+1)-x(:,k+1)) > 1e-5)
                    %     disp('shifting scheme is not correct')
                    %     keyboard
                    % end
                    
                    Px(:,:,k+1) = zeros(obj.nx);
                    for jj = 1:n_sampling_point
                        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    end
                    c = c + kappa*trace(Px(:,:,k+1));
                    
                    % for jj = 1:n_sampling_point
                    % V_comp(:,:,k+1) = (Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,jj,k+1)-x_mean(:,k+1))';
                    % c = c + kappa*trace(w*V_comp(:,:,k+1));
                    % V_covariance(:,:,k+1) = V_covariance(:,:,k+1) + w*V_comp(:,:,k+1);
                    
                    % % debugging
                    % c_variance(j,k+1) = kappa*trace(w*V_comp);
                    %
                    % V_comp_x = (Sig(1,j,k+1)-x_mean(1,k+1))*(Sig(1,j,k+1)-x_mean(1,k+1))';
                    % c_variance_x(j,k+1) = kappa*trace(w*V_comp_x);
                    % V_comp_xd = (Sig(7,j,k+1)-x_mean(7,k+1))*(Sig(7,j,k+1)-x_mean(7,k+1))';
                    % c_variance_xd(j,k+1) = kappa*trace(w*V_comp_xd);
                    %
                    % V_comp_z = (Sig(3,j,k+1)-x_mean(3,k+1))*(Sig(3,j,k+1)-x_mean(3,k+1))';
                    % c_variance_z(j,k+1) = kappa*trace(w*V_comp_z);
                    % V_comp_zd = (Sig(9,j,k+1)-x_mean(9,k+1))*(Sig(9,j,k+1)-x_mean(9,k+1))';
                    % c_variance_zd(j,k+1) = kappa*trace(w*V_comp_zd);
                    % end
                    
                    % accumulate returned cost
                    %c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;%i.i.d mean deviation version
                    %c = c + 1/2*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)) + 1/2*log(det(Px(:,:,k+1)));%ML mean deviation version
                    % % debugging
                    % c_quadratic(k+1) = norm(x(:,k+1)-x_mean(:,k+1))^2;
                    % c_quadratic_x(k+1) = norm(x(1,k+1)-x_mean(1,k+1))^2;
                    % c_quadratic_xd(k+1) = norm(x(7,k+1)-x_mean(7,k+1))^2;
                    % c_quadratic_z(k+1) = norm(x(3,k+1)-x_mean(3,k+1))^2;
                    % c_quadratic_zd(k+1) = norm(x(7,k+1)-x_mean(9,k+1))^2;
                    
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
                        
                        tic
                        % gradient w.r.t state x
                        dSig_m_kplus1_dx_sum = zeros(obj.nx);
                        % gradient w.r.t control u
                        dSig_m_kplus1_du_sum = zeros(obj.nx,nu);
                        dSig_i_kplus1_dx_resample = zeros(obj.nx,obj.nx,n_sampling_point);
                        dSig_i_kplus1_du_resample = zeros(obj.nx,nu,n_sampling_point);
                        
                        for i=1:n_sampling_point
                            if i == 1
                                for m = 1:(2*(obj.nx+nw))% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                                    % gradient of Tr(V_{k+1}) w.r.t control x and u
                                    dSig_m_kplus1_dx = zeros(obj.nx);
                                    dSig_m_kplus1_du = zeros(obj.nx,1);
                                    
                                    chain_rule_indx = k-j;
                                    if j ~= 1
                                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,j+1);
                                    else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2);
                                    end
                                    dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                                    
                                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                        dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                                        dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;

                                        chain_rule_indx = chain_rule_indx - 1;
                                    end
%                                     while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
%                                         dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_dx;
%                                         dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*eta_final(k+2-chain_rule_indx)*(1-w_avg)*dSig_m_kplus1_du;
%                                         chain_rule_indx = chain_rule_indx - 1;
%                                     end
                                    dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                                    dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                                end
                            end
                            
                            % run 2*(obj.nx+nw) times in total to obtain
                            % gradient w.r.t sigma points
                            dSig_i_kplus1_dx = zeros(obj.nx);
                            dSig_i_kplus1_du = zeros(obj.nx,nu);
                            chain_rule_indx = k-j;
                            if j ~= 1
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1);% + dfdSig(:,:,i,j+1);
                            else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                                dSig_i_kplus1_dx = dfdx(:,:,i,j+1) + dfdSig(:,:,i,2)*eye(obj.nx);
                            end
                            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
                            
                            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                                chain_rule_indx = chain_rule_indx - 1;
                            end
                            
                            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
                            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                            
%                             % new sigma point due to resampling mechanism
%                             dSig_i_kplus1_dx_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_dx - w_avg*dSig_m_kplus1_dx_sum);
%                             dSig_i_kplus1_du_resample(:,:,i) = eta_final(k+1)*(dSig_i_kplus1_du - w_avg*dSig_m_kplus1_du_sum);
                        end
                        
%                         dSig_kplus1_dx_sum_resample = zeros(obj.nx,obj.nx);
%                         dSig_kplus1_du_sum_resample = zeros(obj.nx,nu);
%                         
%                         for i =1:n_sampling_point
%                             dSig_kplus1_dx_sum_resample = dSig_kplus1_dx_sum_resample + dSig_i_kplus1_dx_resample(:,:,i);
%                             dSig_kplus1_du_sum_resample = dSig_kplus1_du_sum_resample + dSig_i_kplus1_du_resample(:,:,i);
%                         end
%                         
%                         for i =1:n_sampling_point
%                             dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx_resample(:,:,i) - w_avg*dSig_kplus1_dx_sum_resample);
%                             dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du_resample(:,:,i) - w_avg*dSig_kplus1_du_sum_resample);
%                         end
                        
                        % gradient of mean residual w.r.t state x and control u, assume norm 2
                        dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_dx_sum);%i.i.d mean deviation version
                        dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_avg*dSig_m_kplus1_du_sum);%i.i.d mean deviation version
                        %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_dx_sum);%ML mean deviation version
                        %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + (pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1)))'*(-w_avg*dSig_m_kplus1_du_sum);%ML mean deviation version
                        
                        % dCovdx = zeros(nx,nx,nx);
                        % dCovdu = zeros(nx,nx,nu);
                        % toc
                        %
                        % tic
                        % for pp=1:nx
                        %     for i=1:n_sampling_point
                        %         dCovdmeandev = zeros(nx,obj.nx,nx);%d covariance /d (mean deviation)
                        %         for nn=1:nx
                        %             dCovdmeandev(nn,nn,nn) = 2*w*(Sig(nn,i,k+1)-x_mean(nn,k+1));
                        %             if nn == 1
                        %                 dCovdmeandev(1,2:end,nn) = w*(Sig(2:nx,i,k+1)-x_mean(2:end,k+1))';
                        %                 dCovdmeandev(2:end,1,nn) = dCovdmeandev(1,2:end,nn)';
                        %             elseif nn == nx
                        %                 dCovdmeandev(nx,1:end-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                        %                 dCovdmeandev(1:end-1,nx,nn) = dCovdmeandev(nx,1:end-1,nn)';
                        %             else
                        %                 dCovdmeandev(nn,1:nn-1,nn) = w*(Sig(1:nn-1,i,k+1)-x_mean(1:nn-1,k+1))';
                        %                 dCovdmeandev(nn,nn+1:end,nn) = w*(Sig(nn+1:nx,i,k+1)-x_mean(nn+1:end,k+1))';
                        %                 dCovdmeandev(:,nn,nn) = dCovdmeandev(nn,:,nn)';
                        %             end
                        %             dCovdx(:,:,pp) = dCovdx(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_dx_resample(nn,pp,i) - w_avg*dSig_kplus1_dx_sum_resample(nn,pp));
                        %             if pp <= nu
                        %                 dCovdu(:,:,pp) = dCovdu(:,:,pp) + dCovdmeandev(:,:,nn)*(dSig_i_kplus1_du_resample(nn,pp,i) - w_avg*dSig_kplus1_du_sum_resample(nn,pp));
                        %             end
                        %         end
                        %     end
                        %     dmeanRdx(j,pp,k+1) = dmeanRdx(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdx(:,:,pp));
                        %     if pp <= nu
                        %         dmeanRdu(j,pp,k+1) = dmeanRdu(j,pp,k+1) - 1/2*trace(pinv(Px(:,:,k+1))*(x(:,k+1)-x_mean(:,k+1))*(x(:,k+1)-x_mean(:,k+1))'*pinv(Px(:,:,k+1))*dCovdu(:,:,pp));
                        %     end
                        % end
                        % toc
                    end
                end
                
                dc = [];
                % cost gradient w.r.t x at first time step is zero
                for jj=1:obj.N % index for x_k
                    dTrV_sum_dx_k = zeros(1, obj.nx);
                    if (jj == 1)
                        dmeanR_sum_dx_k = 2*(x(:,jj)-x_mean(:,jj))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
                    else
                        dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';%i.i.d mean deviation version
                        %dmeanR_sum_dx_k = (inv(Px(:,:,jj))*(x(:,jj)-x_mean(:,jj)))';%ML mean deviation version
                    end
                    for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(jj,:,jjj);
                        dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(jj,:,jjj);
                    end
                    %dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
                    %dc = [dc, dmeanR_sum_dx_k];
                    dc = [dc, kappa*dTrV_sum_dx_k];
                    
                end
                
                % cost gradient w.r.t u at first time step is zero, since
                for jj=1:obj.N % index for u_k
                    dTrV_sum_du_k = zeros(1, nu);
                    dmeanR_sum_du_k = zeros(1, nu);
                    for jjj = jj+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
                        dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(jj,:,jjj);
                        dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(jj,:,jjj);
                    end
                    %dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
                    %dc = [dc, dmeanR_sum_du_k];
                    dc = [dc, kappa*dTrV_sum_du_k];
                end
                
                % scale this robust cost
                c = obj.options.contact_robust_cost_coeff*c;
                dc = obj.options.contact_robust_cost_coeff*dc;
            end
        end
        
        function [c,dc] = getTimeStep(obj, h)
            global timestep_updated
            timestep_updated = h;
            c = 0;
            dc = 0;
        end
        
        function [c,dc] = robustLCPcost(obj, slack_var)
            c = obj.options.robustLCPcost_coeff*sum(slack_var);
            dc = obj.options.robustLCPcost_coeff*ones(1,length(slack_var));
            fprintf('sum of slack variable cost: %4.4f\n',c);
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
            
            BuminusC = B*u-C + [Fext(1);0;Fext(2);zeros(3,1)];% manually add an external force
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
                    dfqdFext = zeros(nq,2);
                    dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl), dfqdFext];% expand the dimension with Fext
            end
            
            % H*v1 = H*v0 + h*(B*u - C) + n^T lambda_N + d^T * lambda_f
            fv = H*(v1 - v0) - h*BuminusC;
            % [h q0 v0 q1 v1 u l ljl]
            
            %dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl)] + ...
            %    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,nu+nl+njl)];
            % expand the dimension with Fext
            dfvdFext = -h*[1,0;zeros(1,2);0,1;zeros(3,2)];
            dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl), dfvdFext] + ...
                [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,2+nu+nl+njl)];
            
            if nl>0
                global uncertain_mu;
                global terrain_index;
                manip = obj.plant.getManipulator();
                if strcmp(manip.uncertainty_source, 'friction_coeff')
                    uncertain_mu = [];%reset mu, and take the average
                    [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                elseif strcmp(manip.uncertainty_source, 'terrain_height')
                    sample_num = length(manip.uncertain_position_set);
                    for ii=1:sample_num
                        manip.uncertain_phi = manip.uncertain_position_set(:,ii);
                        terrain_index = ii;
                        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q1,false,obj.options.active_collision_options);
                        phi_stack(:,ii) = phi;
                        n_stack(:,:,ii) = n;
                        dn_stack(:,:,ii) = dn;
                        D_stack1(:,:,ii) = D{1};
                        D_stack2(:,:,ii) = D{2};
                        D_stack3(:,:,ii) = D{3};
                        D_stack4(:,:,ii) = D{4};
                        dD_stack1(:,:,ii) = dD{1};
                        dD_stack2(:,:,ii) = dD{2};
                        dD_stack3(:,:,ii) = dD{3};
                        dD_stack4(:,:,ii) = dD{4};
                    end
                    phi = sum(phi_stack,2)/sample_num;
                    n = sum(n_stack,3)/sample_num;
                    dn = sum(dn_stack,3)/sample_num;
                    D{1} = sum(D_stack1,3)/sample_num;
                    D{2} = sum(D_stack2,3)/sample_num;
                    D{3} = sum(D_stack3,3)/sample_num;
                    D{4} = sum(D_stack4,3)/sample_num;
                    dD{1} = sum(dD_stack1,3)/sample_num;
                    dD{2} = sum(dD_stack2,3)/sample_num;
                    dD{3} = sum(dD_stack3,3)/sample_num;
                    dD{4} = sum(dD_stack4,3)/sample_num;
                elseif strcmp(manip.uncertainty_source, 'friction_coeff+terrain_height')
                    sample_num = length(manip.uncertain_mu_set);
                    for ii=1:sample_num
                        manip.uncertain_mu = manip.uncertain_mu_set(ii);% uncertainy in mu is streamed down into the contactConstraint_manual function
                        manip.uncertain_phi = manip.uncertain_position_set(:,ii);
                        uncertain_mu = manip.uncertain_mu_set(ii);
                        terrain_index = ii;
                        [phi,~,~,~,~,~,~,~,n,D,dn,dD] = manip.plant_sample{terrain_index}.contactConstraints(q1,false,obj.options.active_collision_options);
                        phi_stack(:,ii) = phi;
                        n_stack(:,:,ii) = n;
                        dn_stack(:,:,ii) = dn;
                        D_stack1(:,:,ii) = D{1};
                        D_stack2(:,:,ii) = D{2};
                        D_stack3(:,:,ii) = D{3};
                        D_stack4(:,:,ii) = D{4};
                        dD_stack1(:,:,ii) = dD{1};
                        dD_stack2(:,:,ii) = dD{2};
                        dD_stack3(:,:,ii) = dD{3};
                        dD_stack4(:,:,ii) = dD{4};
                    end
                    phi = sum(phi_stack,2)/sample_num;
                    n = sum(n_stack,3)/sample_num;
                    dn = sum(dn_stack,3)/sample_num;
                    D{1} = sum(D_stack1,3)/sample_num;
                    D{2} = sum(D_stack2,3)/sample_num;
                    D{3} = sum(D_stack3,3)/sample_num;
                    D{4} = sum(D_stack4,3)/sample_num;
                    dD{1} = sum(dD_stack1,3)/sample_num;
                    dD{2} = sum(dD_stack2,3)/sample_num;
                    dD{3} = sum(dD_stack3,3)/sample_num;
                    dD{4} = sum(dD_stack4,3)/sample_num;
                elseif isempty(manip.uncertainty_source)%non-robust version
                    [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                end
                %[phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
                
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
                
                fv = fv - J'*lambda;% lambda is the impulse, i.e., force*h, h is time step.
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
            
            % f_numeric = f;
            % df_numeric = df;
            % disp('check gradient')
            % [f_numeric,df_numeric] = geval(@(h,x0,x1,lambda,Fext) dynamics_constraint_fun_check(obj,h,x0,x1,u,lambda,lambda_jl,Fext),h,x0,x1,lambda,Fext,struct('grad_method','numerical'));
            % valuecheck(df,df_numeric,1e-3);
            % valuecheck(f,f_numeric,1e-3);
            
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
                
                BuminusC = B*u-C + [Fext(1);0;Fext(2);zeros(3,1)];% manually add an external force
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
                        dfqdFext = zeros(nq,2);
                        dfq = [-v1, -eye(nq), zeros(nq,nv), eye(nq), -h*eye(nq), zeros(nq,nu+nl+njl), dfqdFext];% expand the dimension with Fext
                end
                
                % H*v1 = H*v0 + h*(B*u - C) + n^T lambda_N + d^T * lambda_f
                fv = H*(v1 - v0) - h*BuminusC;
                % [h q0 v0 q1 v1 u l ljl]
                
                %dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl)] + ...
                %    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,nu+nl+njl)];
                % expand the dimension with Fext
                dfvdFext = -h*[1,0;zeros(1,2);0,1;zeros(3,2)];
                dfv = [-BuminusC, zeros(nv,nq), -H, zeros(nv,nq), H,-h*B, zeros(nv,nl+njl), dfvdFext] + ...
                    [zeros(nv,1) matGradMult(dH0,v1-v0)-h*dBuminusC0 matGradMult(dH1,v1-v0)-h*dBuminusC1 zeros(nv,2+nu+nl+njl)];
                
                if nl>0
                    [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(q1,false,obj.options.active_collision_options);
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
            Fexttraj = PPTrajectory(foh(t,[reshape(z(obj.F_ext_inds),[], obj.N)]));
        end
        
        function obj = setupVariables(obj,N)
            obj = setupVariables@DirectTrajectoryOptimization_Brick(obj,N);
            [~,normal,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), false, obj.options.active_collision_options);
            obj.nC = size(normal, 2);
            obj.nD = 2*length(d);
            obj.nFext = 2;%[set up for brick falling example]
            
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
            
            obj.F_ext_inds = reshape(obj.num_vars + (1:N*obj.nFext), obj.nFext, N);
            obj = obj.addDecisionVariable(N*obj.nFext);
            
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
                z0(obj.F_ext_inds) = F_ext_eval(:,1:obj.N); % take N*obj.nFext values
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
    end
end