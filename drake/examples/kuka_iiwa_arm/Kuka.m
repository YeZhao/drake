classdef Kuka < RigidBodyManipulator
    properties
        curr_state
        goal_state
        dircol_planner
    end
    
    methods
        function obj = Kuka()
            % suppress the cylinder warning
            w = warning('off','all');
            
            % load the object
            urdf_path = [getDrakePath,'/examples/kuka_iiwa_arm/urdf/iiwa14.urdf'];
            options.floating = false;
            obj@RigidBodyManipulator(urdf_path, options);
            
            warning(w);
        end
    end
    
    methods(Static = true)
        function [h,dh] = finalCost(t,x,goal_state)
            % simple quadratic cost to the goal
            Q = 10*eye(length(x));
            T = 1;

            err = x - goal_state;
            t_err = t;
            h = T*t_err + err'*Q*err;
            dh = [T,2*err'*Q];
            % disp('Final Cost')
            %h
         end
        
        function [g,dg] = runningCost(dt,x,u,goal_state)
            % simple quadratic cost to the goal
            R = 0.01*eye(length(u));
            pos_cost = 10*ones(length(x)/2,1);
            vel_cost = 3*ones(length(x)/2,1);
            Q = diag([pos_cost; vel_cost]);

            err = x - goal_state;
            g = err'*Q*err + u'*R*u;
            dg = [0, 2*err'*Q, 2*u'*R];
            % disp('Running Cost')
            %g
        end
    end
end