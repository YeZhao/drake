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
        
        function handle_request(obj, msg)
            disp('Handling plan request')
            msg = drake.lcmt_matlab_plan_request(msg);
            [t,x,u,info] = obj.plan_dircol(msg.start_state, msg.goal_state);
            res = drake.lcmt_matlab_plan_response();
            res.num_timesteps = length(t);
            res.state_size = 14;
            res.input_size = 7;
            res.time = t;
            disp('start state')
            msg.start_state
            disp('end_state')
            msg.goal_state
            
            res.state = x;
            res.input = u;
            res.info = info;
  
           obj.lc.publish(obj.response_channel, res) 
        end
        
        function [t, x, u, info] = plan_dircol(obj, x0, xG)
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
            prog = DircolTrajectoryOptimization(obj.kuka, N, [2 6]);% arbitrary timescale
            prog.addStateConstraint(ConstantConstraint(x0),1);
            prog.addStateConstraint(ConstantConstraint(xG),N);
            prog.addRunningCost(@obj.kuka.runningCost);
            prog.addFinalCost(@obj.kuka.finalCost);

            % initial guess at the trajectory
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xG)]));

            % solve with dircol
            [xtraj, utraj, z, F, info] = prog.solveTraj(tf0,traj_init);

            % unpack the trajectory
            t = linspace(xtraj.tspan(1),xtraj.tspan(end),N);
            x = zeros(14,N);
            u = zeros(7,N);
            for i=1:N
                x(:,i) = xtraj.eval(t(i));
                u(:,i) = utraj.eval(t(i));
            end
        end
        
        function run(obj)
            disp('In the main loop')
            while true
                % check if there is a message available
                req_msg = obj.lcmAggregator.getNextMessage(5);
                if isempty(req_msg)
                    continue
                end
                obj.handle_request(req_msg);
            end
        end
    end
end