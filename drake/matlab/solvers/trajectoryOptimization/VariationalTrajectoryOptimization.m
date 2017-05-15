classdef VariationalTrajectoryOptimization < DirectTrajectoryOptimization
    
    properties
        nC %Number of contact points
        nD %Number of basis vectors in polyhedral friction cone
        nL %Number of contact constraints
        twoD %Is the plant a PlanarRigidBodyManipulator?
        IntegrationMethod %Midpoint or Simpson's rule for integration
        l_inds; %Indices of all contact-related variables
    end
    
    properties (Constant)
        MIDPOINT = 2;
        SIMPSON = 3;
    end
    
    methods
        function obj = VariationalTrajectoryOptimization(plant,N,duration,options)
            if nargin<4, options=struct(); end
            
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
            end
            
            typecheck(plant,'RigidBodyManipulator');
            if isa(plant,'PlanarRigidBodyManipulator')
                obj.twoD = true;
            end
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
        end
        
        function obj = setupVariables(obj, N)
            switch obj.options.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.MIDPOINT;
                    
                    [phi,~,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), false, obj.options.active_collision_options);
                    
                    Nh = N-1;
                    Nq = obj.plant.getNumPositions();
                    nU = obj.plant.getNumInputs();
                    nC = length(phi);
                    nD = 2*length(d);
                    nL = (3+2*nD)*nC; %Includes contact forces, contact distance, tangential speed, and friction cone
                    
                    obj.N = N;
                    obj.nC = nC;
                    obj.nD = nD;
                    obj.nL = nL;
                    
                    num_vars = Nh + N*Nq + (N-1)*nU + (N-1)*nL;
                    
                    obj.h_inds = (1:Nh)';
                    obj.x_inds = reshape(Nh + (1:(N*Nq)), Nq, N);
                    obj.u_inds = reshape(Nh+N*Nq + (1:N*nU), nU, N);
                    obj.l_inds = reshape(Nh+N*Nq+(N-1)*nU + (1:((N-1)*nL)), nL, N);
                    
                    x_names = cell(num_vars,1);
                    for i = 1:N
                        if(i<N)
                            x_names{i} = sprintf('h[%d]',i);
                            for j = 1:nU
                                x_names{nH+nX*N+(i-1)*nU+j} = sprintf('u%d[%d]',j,i);
                            end
                            for j = 1:nL
                                x_names{nH+nX*N+nU*(N-1)+(i-1)*nL+j} = sprintf('l%d[%d]',j,i);
                            end
                        end
                        for j = 1:nX
                            x_names{nH+(i-1)*nX+j}=sprintf('x%d[%d]',j,i);
                        end
                    end
                    
                    obj = obj.addDecisionVariable(num_vars,x_names);
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.SIMPSON;
                    
                    
            end
        end
        
        function obj = addDynamicConstraints(obj)
            Nq = obj.plant.getNumPositions();
            Nv = obj.plant.getNumVelocities;
            nU = obj.plant.getNumInputs();
            nC = obj.nC;
            nD = obj.nD;
            N = obj.N;
            assert(Nq == Nv) % not quite ready for the alternative
            
            switch obj.IntegrationMethod
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                case VariationalTrajectoryOptimization.SIMPSON
            
            end
            
        end
        
        function obj = addStateConstraint(obj,constraint,time_index,x_indices)
            
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
                      
                      zeros(nQ, 2*nU+nC+(nC*nD))];
            
            %Contact stuff
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin2 = obj.plant.doKinematics((q2+q3)/2, (q3-q2)/h2, kinopts);
            [~,~,~,~,~,~,~,~,n2,D2,dn2,dD2] = obj.plant.contactConstraints(kin2, obj.multiple_contacts);
            D2 = reshape(cell2mat(D2')',nQ,nC*nD)';
            dD2 = reshape(cell2mat(dD2)',nQ,nC*nD*nD)';
            
            
            %Total dynamics residual incluing control + contact forces
            f = f_del + (h1/2)*B1*u1 + (h2/2)*B2*u2 + h2*(n2*c2 + D2*b2);
            
            %Derivatives
            df = df_del + [0.5*B1*u1, 0.5*B2*u2 + n2*c2 + D2*b2, (h1/2)*kron(u1', eye(nQ))*dB1, (h2/2)*kron(u2', eye(nQ))*dB2, 
            
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
        
        function [xtraj,utraj,ltraj,z,F,info] = solveTraj(obj,t_init,traj_init)
            [xtraj,utraj,z,F,info] = solveTraj@DirectTrajectoryOptimization(obj,t_init,traj_init);
            t = [0; cumsum(z(obj.h_inds))];
            if obj.nC>0
                ltraj = PPTrajectory(foh(t,reshape(z(obj.l_inds),[],obj.N)));
            else
                ltraj = [];
            end
        end
        
        function xtraj = reconstructStateTrajectory(obj,z)
            Nq = obj.plant.getNumPositions();
            
            t = [0; cumsum(z(obj.h_inds))];
            x = reshape(z(obj.x_inds),[],obj.N);
            
            switch obj.IntegrationMethod
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    q = x(1:Nq,:);
                    qtraj = PPTrajectory(foh(tq,q));
                    
                    v = [diff(x,1,2)/z(obj.h_inds(1)), zeros(Nq,1)]; %zoh (correctly in this case) ignores the last entry in v
                    vtraj = PPTrajectory(zoh(t,v));
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
            
            xtraj = [qtraj; vtraj];
            xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
        end
        
        function z0 = getInitialVars(obj,t_init,traj_init)
            
        end
    end
end