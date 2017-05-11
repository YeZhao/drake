classdef VariationalTrajectoryOptimization < DirectTrajectoryOptimization
    
    properties
        Nc %Number of contact points
        Nd %Number of basis vectors in polyhedral friction cone
        IntegrationMethod %Midpoint or Simpson's rule for integration
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
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
        end
        
        function obj = setupVariables(obj, N)
            switch obj.options.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.MIDPOINT;
                    
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.SIMPSON;
                    
                    
            end
        end
        
        function obj = addDynamicConstraints(obj)
            Nq = obj.plant.getNumPositions();
            Nv = obj.plant.getNumVelocities;
            Nu = obj.plant.getNumInputs();
            Nc = obj.Nc;
            Nd = obj.Nd;
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
        
        function [f,df] = midpoint_dynamics_constraint_fun(obj,h,q0,q1,u,lambda)
            Nq = obj.plant.getNumPositions;
            Nu = obj.plant.getNumInputs;
            
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
        
        function z0 = getInitialVars(obj,t_init,traj_init)
            
        end
    end
end