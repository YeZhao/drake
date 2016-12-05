classdef KukaPlanner
    properties
        lc
        kuka
        lcmAggregator
    end
    properties(Constant = true)
        request_channel = 'MATLAB_KUKA_DIRCOL_REQUEST'
        response_channel = 'MATLAB_KUKA_DIRCOL_RESPONSE'
    end
    methods
        function obj = KukaPlanner()
            % create LCM object
            obj.lc = lcm.lcm.LCM.getSingleton();
            
            % subscribe to the MATLAB_KUKA_DIRCON channel
            obj.lcmAggregator = drake.matlab.util.MessageMonitor(drake.lcmt_matlab_plan_request,'timestamp');
            obj.lc.subscribe(obj.request_channel,obj.lcmAggregator)

            % create kuka object and supress warnings
            w = warning('off','all');
            obj.kuka = Kuka;
            warning(w);
        end
        
        function handle_request(obj, request)
            obj.plan_dircol();
        end
        
        function [t, x, u] = plan_dircol(obj, x0, xG)
            % test state and goal
            if (nargin == 0)
                x0 = [1; zeros(6,1); zeros(7,1)];
                xG = zeros(14,1);
            end

            % trajectory params
            tf0 = 4;
            N = 20;
            obj.kuka.goal_state = xG;

            % setup dircol
            prog = DircolTrajectoryOptimization(k, N, [2 6]);% arbitrary timescale
            prog.addStateConstraint(ConstantConstraint(x0),1);
            prog.addStateConstraint(ConstantConstraint(xG),1);
            prog.addRunningCost(@k.runningCost);
            prog.addFinalCost(@k.finalCost);

            % initial guess at the trajectory
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xG)]));

            % solve with dircol
            [xtraj, utraj, z, F, info] = prog.solveTraj(tf0,traj_init);

            % unpack the trajectory
            t = linspace(xtraj.tspan(1),xtraj.tspan(1),N);
            x = zeros(14,N);
            u = zeros(7,N);
            for i=1:N
                x(:,i) = xtraj.eval(t(i));
                u(:,i) = utraj.eval(t(i));
            end
        end
        
        function run(obj)
            while true
                % check if there is a message available
                req_msg = obj.lcmAggregator.getNextMessage(5);
                if isempty(req_msg)
                    continue
                end
                plan = obj.handle_request(req_msg);
                % TODO: implement handle_request method
%                 plan = plan.toLCM();
%                 obj.lc.publish(obj.response_channel, plan)
            end
        end
    end
end