classdef VariationalTrajectoryOptimization < DirectTrajectoryOptimization
    
    properties
        nC %Number of contact points
        nD %Number of basis vectors in polyhedral friction cone
        nL %Number of contact constraints
        integration_method %Midpoint or Simpson's rule for integration
        l_inds; %Indices of all contact-related variables
    end
    
    properties (Constant)
        MIDPOINT = 2;
        SIMPSON = 3;
    end
    
    methods
        function obj = VariationalTrajectoryOptimization(plant,N,duration,options)
            if nargin < 4
                options=struct(); 
            end
            
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
            end
            if ~isfield(options,'multiple_contacts')
                options.multiple_contacts = false;
            end
            
            typecheck(plant,'RigidBodyManipulator');
            if isa(plant,'PlanarRigidBodyManipulator')
                options.twoD = true;
            else
                options.twoD = false;
            end
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);

        end
        
        function obj = setupVariables(obj, N)
            switch obj.options.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    obj.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
                    
                    [phi,~,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), false, obj.options.active_collision_options);
                    
                    nH = N-1;
                    nQ = obj.plant.getNumPositions();
                    nU = obj.plant.getNumInputs();
                    nC = length(phi);
                    nD = 2*length(d);
                    nL = 1+(3+2*nD)*nC; %Includes contact forces, contact distance, tangential speed, friction cone, and smoothing parameter
                    
                    obj.N = N;
                    obj.nC = nC;
                    obj.nD = nD;
                    obj.nL = nL;
                    
                    num_vars = nH + N*nQ + (N-1)*nU + (N-2)*nL;
                    
                    obj.h_inds = (1:nH)';
                    obj.x_inds = reshape(nH + (1:(N*nQ)), nQ, N);
                    obj.u_inds = reshape(nH+N*nQ + (1:((N-1)*nU)), nU, N-1);
                    obj.l_inds = reshape(nH+N*nQ+(N-1)*nU + (1:((N-2)*nL)), nL, N-2);
                    
                    x_names = cell(num_vars,1);
                    for i = 1:N
                        if(i<N)
                            x_names{i} = sprintf('h[%d]',i);
                            for j = 1:nU
                                x_names{nH+nQ*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
                            end
                        end
                        if(i<(N-1))
                            for j = 1:nL
                                x_names{nH+nQ*N+nU*(N-1)+(i-1)*nL+j} = sprintf('l%d[%d]',j,i);
                            end
                        end
                        for j = 1:nQ
                            x_names{nH+(i-1)*nQ+j}=sprintf('q%d[%d]',j,i);
                        end
                    end
                    
                    obj = obj.addDecisionVariable(num_vars,x_names);
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.SIMPSON;
                    
                    
            end
        end
        
        function obj = addDynamicConstraints(obj)
            nQ = obj.plant.getNumPositions();
            nV = obj.plant.getNumVelocities();
            nU = obj.plant.getNumInputs();
            nC = obj.nC;
            nD = obj.nD;
            N = obj.N;
            assert(nQ == nV) % can't handle quaternions yet
            
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    dyn_constraints = cell(N-2,1);
                    dyn_inds = cell(N-2,1);
                    
                    cont_constraints = cell(N-2,1);
                    cont_inds = cell(N-2,1);
                    
                    ineq_constraints = cell(N-2,1);
                    ineq_inds = cell(N-2,1);
                    
                    nvars1 = 2+3*nQ+2*nU+nC+nD*nC;
                    cnstr1 = FunctionHandleConstraint(zeros(nQ,1), zeros(nQ,1), nvars1, @obj.midpoint_dynamics_constraint_fun);
                    
                    nvars2 = 1+2*nQ+2*nC+nD*nC+nC+nD*nC+1;
                    cnstr2 = FunctionHandleConstraint(zeros(3*nC+3*nD*nC+1,1), zeros(3*nC+3*nD*nC+1,1), nvars2, @obj.midpoint_contact_constraint_fun);
                    
                    nvars3 = obj.nL;
                    cnstr3 = BoundingBoxConstraint(zeros(obj.nL,1), inf*ones(obj.nL,1));
                    
                    for i = 1:obj.N-2
                        dyn_inds{i} = {obj.h_inds(i); obj.h_inds(i+1); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.x_inds(:,i+2); obj.u_inds(:,i); obj.u_inds(:,i+1); obj.l_inds(2*nC+nD*nC+(1:nC),i); obj.l_inds(3*nC+nD*nC+(1:nD*nC),i)};
                        dyn_constraints{i} = cnstr1;
                        obj = obj.addConstraint(dyn_constraints{i}, dyn_inds{i});
                        
                        cont_inds{i} = {obj.h_inds(i+1); obj.x_inds(:,i+1); obj.x_inds(:,i+2); obj.l_inds(1:nC,i); obj.l_inds(nC+(1:nC),i); obj.l_inds(2*nC+(1:nD*nC),i); obj.l_inds(2*nC+nD*nC+(1:nC),i); obj.l_inds(3*nC+nD*nC+(1:nD*nC),i); obj.l_inds(end,i)};
                        cont_constraints{i} = cnstr2;
                        obj = obj.addConstraint(cont_constraints{i}, cont_inds{i});
                        
                        ineq_inds{i} = {obj.l_inds(:,i)};
                        ineq_constraints{i} = cnstr3;
                        obj = obj.addConstraint(ineq_constraints{i}, ineq_inds{i});
                    end
                    
                    
                case VariationalTrajectoryOptimization.SIMPSON
            
            end
            
        end
        
        function obj = addStateConstraint(obj,constraint,time_index,x_indices)
            if ~iscell(time_index)
                % then use { time_index(1), time_index(2), ... } ,
                % aka independent constraints for each time
                time_index = num2cell(reshape(time_index,1,[]));
            end
            if nargin < 4
                x_indices = 1:size(obj.x_inds,1);
            end
            
            for j=1:length(time_index),
                cstr_inds = mat2cell(obj.x_inds(x_indices,time_index{j}),numel(x_indices),ones(1,length(time_index{j})));
                
                % record constraint for posterity
                obj.constraints{end+1}.constraint = constraint;
                obj.constraints{end}.var_inds = cstr_inds;
                obj.constraints{end}.time_index = time_index;
                
                obj = obj.addConstraint(constraint,cstr_inds);
            end
        end
        
        function obj = addRunningCost(obj,running_cost_function)
            
        end
        
        function obj = addInitialCost(obj,initial_cost)
            
        end
        
        function obj = addFinalCost(obj,final_cost_function)
            
        end
        
        function [f,df] = midpoint_dynamics_constraint_fun(obj,h1,h2,q1,q2,q3,u1,u2,c2,b2)
            
            nC = obj.nC;
            nD = obj.nD;
            nQ = length(q1);
            nU = length(u1);
            
            %Discrete Euler-Lagrange equation
            [D1L1,D2L1,D1D1L1,D1D2L1,D2D2L1,B1,dB1] = obj.LagrangianDerivs((q1+q2)/2,(q2-q1)/h1);
            [D1L2,D2L2,D1D1L2,D1D2L2,D2D2L2,B2,dB2] = obj.LagrangianDerivs((q2+q3)/2,(q3-q2)/h2);
            f_del = (h1/2)*D1L1 + D2L1 + (h2/2)*D1L2 - D2L2;
            
            df_del = [0.5*D1L1, 0.5*D1L2, ... % d/dh1, d/dh2
                      (h1/2)*((1/2)*D1D1L1-(1/h1)*D1D2L1)+(1/2)*D1D2L1-(1/h1)*D2D2L1, ... % d/dq1
                      (h1/2)*((1/2)*D1D1L1+(1/h1)*D1D2L1)+(1/2)*D1D2L1+(1/h1)*D2D2L1 + (h2/2)*((1/2)*D1D1L2-(1/h2)*D1D2L2)-(1/2)*D1D2L2+(1/h2)*D2D2L2, ... % d/dq2
                      (h2/2)*((1/2)*D1D1L2+(1/h2)*D1D2L2) - (1/2)*D1D2L2 - (1/h2)*D2D2L2, ... % d/dq3
                      zeros(nQ, 2*nU+nC+(nC*nD))];
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin2 = obj.plant.doKinematics((q2+q3)/2, (q3-q2)/h2, kinopts);
            [~,~,~,~,~,~,~,~,n2,D2,dn2,dD2] = obj.plant.contactConstraints(kin2, obj.options.multiple_contacts);
            D2 = reshape(cell2mat(D2')',nQ,nC*nD)';
            dD2 = reshape(cell2mat(dD2)',nQ,nC*nD*nQ)';
            
            %Total dynamics residual incluing control + contact forces
            f = f_del + (h1/2)*B1*u1 + (h2/2)*B2*u2 + h2*(n2'*c2 + D2'*b2);
            
            %Dynamics Derivatives
            df = df_del + [(1/2)*B1*u1, (1/2)*B2*u2 + n2'*c2 + D2'*b2, ... % d/dh1, d/dh2
                           (h1/4)*kron(u1', eye(nQ))*dB1, ... % d/dq1
                           (h1/4)*kron(u1', eye(nQ))*dB1 + (h2/4)*kron(u2', eye(nQ))*dB2 + (h2/2)*kron(c2', eye(nQ))*dn2 + (h2/2)*kron(b2', eye(nQ))*dD2, ... % d/dq2
                           (h2/4)*kron(u2', eye(nQ))*dB2 + (h2/2)*kron(c2', eye(nQ))*dn2 + (h2/2)*kron(b2', eye(nQ))*dD2, ... % d/dq3
                           (h1/2)*B1, (h2/2)*B2, h2*n2', h2*D2']; % d/du1, d/du2, d/dc2, d/db2
        end
        
        function [f,df] = midpoint_contact_constraint_fun(obj,h,q1,q2,phi,psi,eta,c,b,s)
            mu = 1; %This is currently hard-coded in Drake
            nC = obj.nC;
            nD = obj.nD;
            nQ = length(q1);
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.plant.doKinematics(q2, (q2-q1)/h, kinopts);
            [gap,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts);
            D = reshape(cell2mat(D')',nQ,nC*nD)';
            dD = reshape(cell2mat(dD)',nQ,nC*nD*nQ)';
            
            e = ones(nD,1);
            E = kron(eye(nC),e');
            
            %Normal force
            f1 = gap - phi;
            [f2,dfa2,dfb2,dfs2] = obj.smoothfb(phi,c,s);
            
            %Tangential velocity
            [f3,dfa3,dfb3,dfs3] = obj.smoothfb((mu*c - E*b),psi,s);
            
            %Polyhedral friction cone
            f4 = E'*psi + D*((q2-q1)/h) - eta;
            [f5,dfa5,dfb5,dfs5] = obj.smoothfb(E'*phi,b,s);
            [f6,dfa6,dfb6,dfs6] = obj.smoothfb(eta,b,s);
            
            f = [f1; f2; f3; f4; f5; f6; exp(s)-1];
            
            df = [zeros(nC,1), zeros(nC,nQ), n, -eye(nC), zeros(nC,nC), zeros(nC,nD*nC), zeros(nC,nC), zeros(nC,nD*nC), zeros(nC,1);
                  zeros(nC,1), zeros(nC,nQ), zeros(nC,nQ), dfa2, zeros(nC,nC), zeros(nC,nD*nC), dfb2, zeros(nC,nD*nC), dfs2;
                  zeros(nC,1), zeros(nC,nQ), zeros(nC,nQ), zeros(nC,nC), dfb3, zeros(nC,nD*nC), mu*dfa3, -dfa3*E, dfs3;
                  -D*(q2-q1)/(h*h), -D/h, D/h + kron((q2-q1)'/h, eye(nD*nC))*dD, zeros(nD*nC,nC), E', -eye(nD*nC), zeros(nD*nC,nC), zeros(nD*nC,nD*nC), zeros(nD*nC,1);
                  zeros(nD*nC,1), zeros(nD*nC,nQ), zeros(nD*nC,nQ), dfa5*E', zeros(nD*nC,nC), zeros(nD*nC,nD*nC), zeros(nD*nC,nC), dfb5, dfs5;
                  zeros(nD*nC,1), zeros(nD*nC,nQ), zeros(nD*nC,nQ), zeros(nD*nC,nC), zeros(nD*nC,nC), dfa6, zeros(nD*nC,nC), dfb6, dfs6;
                  zeros(1,1), zeros(1,nQ), zeros(1,nQ), zeros(1,nC), zeros(1,nC), zeros(1,nD*nC), zeros(1,nC), zeros(1,nD*nC), exp(s)];
            
        end
        
        function [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dBdq] = LagrangianDerivs(obj,q,v)
            nq = length(q);
            nv = length(v);
            [M,G,B,dM,dG,dB] = manipulatorDynamics(obj.plant, q, zeros(nv,1));
            dM = reshape(dM,nq*nq,nq+nv);
            dMdq = dM(:,1:nq);
            dBdq = dB(:,1:nq);
            D1L = 0.5*dMdq'*kron(v,v) - G;
            D2L = M*v;
            D1D1L = -dG(:,1:nq); %throwing out second derivative of M terms here
            D1D2L = kron(v',eye(nq))*dMdq;
            D2D2L = M;
        end
        
        function [f, dfda, dfdb, dfds] = smoothfb(obj,a,b,s)
            Na = length(a);
            f0 = sqrt(a.*a + b.*b  + s*s*ones(Na,1));
            f = f0 - (a + b);
            dfda = diag(a./f0) - eye(Na);
            dfdb = diag(b./f0) - eye(Na);
            dfds = s*ones(Na,1)./f0;
        end
        
        function [xtraj,utraj,ltraj,z,F,info] = solveTraj(obj,t_init,traj_init)
            [xtraj,utraj,z,F,info] = solveTraj@DirectTrajectoryOptimization(obj,t_init,traj_init);
            t = [0; cumsum(z(obj.h_inds))];
            if obj.nC>0
                ltraj = PPTrajectory(zoh(t(3:end),reshape(z(obj.l_inds),[],obj.N-2)));
            else
                ltraj = [];
            end
        end
        
        function xtraj = reconstructStateTrajectory(obj,z)
            nQ = obj.plant.getNumPositions();
            
            t = [0; cumsum(z(obj.h_inds))];
            x = reshape(z(obj.x_inds),[],obj.N);
            
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    q = x(1:nQ,:);
                    qtraj = PPTrajectory(foh(t,q));
                    
                    v = [diff(x,1,2)/z(obj.h_inds(1)), zeros(nQ,1)]; %zoh (correctly in this case) ignores the last entry in v
                    vtraj = PPTrajectory(zoh(t,v));
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
            
            xtraj = [qtraj; vtraj];
            xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
        end
        
        function z0 = getInitialVars(obj,t_init,traj_init)
            nQ = obj.plant.getNumPositions();
            
            if isscalar(t_init)
                t_init = linspace(0,t_init,obj.N);
            elseif length(t_init) ~= obj.N
                error('The initial sample times must have the same length as property N')
            end
            z0 = zeros(obj.num_vars,1);
            z0(obj.h_inds) = diff(t_init);
            
            if nargin < 3
                traj_init = struct();
            end
            
            nU = getNumInputs(obj.plant);
            if isfield(traj_init,'u')
                z0(obj.u_inds) = traj_init.u.eval(t_init(1:end-1));
            else
                z0(obj.u_inds) = 0.01*randn(nU,obj.N-1);
            end
            
            if isfield(traj_init,'x')
                xsamp = traj_init.x.eval(t_init);
                z0(obj.x_inds) = xsamp(1:nQ,:);
            end
            
            %Set smoothing parameters to 1
            z0(obj.l_inds(end,:)) = 1;
        end
    end
end