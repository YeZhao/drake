classdef Kuka < RigidBodyManipulator
    properties
        curr_state
        goal_state
        dircol_planner
    end
    
    methods
        function obj = Kuka()
            % suppress the cylinder warning
            warning('off','Drake:RigidBodyManipulator:ReplacedCylinder');
            
            % load the object
            urdf_path = [getDrakePath,'/examples/kuka_iiwa_arm/urdf/iiwa14.urdf'];
            options.floating = false;
            obj@RigidBodyManipulator(urdf_path, options);
        end
        
        function [g,dg] = runningCost(dt,x,u)
            % simple quadratic cost to the goal
            err = x - obj.goal_state;
            g = err'*err + u'*u;
            dg = [0,2*err',2*u'];
        end
        
        function [h,dh] = finalCost(t,x)
            % simple quadratic cost to the goal
            err = x - obj.goal_state;
            g = t + err'*err;
            dg = [1,2*err'];
        end
    end
end