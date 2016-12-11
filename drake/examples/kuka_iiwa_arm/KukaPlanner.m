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
            obj.kuka = Kuka;
            
            % subscribe to the MATLAB_KUKA_DIRCON channel
            obj.lcmAggregator = drake.matlab.util.MessageMonitor(drake.lcmt_matlab_plan_request,'timestamp');
            obj.lc.subscribe(obj.request_channel,obj.lcmAggregator)
        end
        
        function handle_request(obj, msg)
            disp('Handling plan request')
            msg = drake.lcmt_matlab_plan_request(msg);
            [t,x,u,info] = obj.plan_dircol(msg.start_state, msg.goal_state);
            res = drake.lcmt_matlab_plan_response();
            res.num_timesteps = length(t);
            res.state_size = obj.kuka.getNumStates;
            res.input_size = obj.kuka.getNumPositions;
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
            Nx = obj.kuka.getNumStates;
            Nq = obj.kuka.getNumPositions;
            Nv = Nx-Nq;
            
            % test state and goal
            if (nargin == 0)
                x0 = [1; zeros(Nq-1,1); zeros(Nv,1)];
                xG = zeros(Nx,1);
            end

            % trajectory params
            tf0 = 2;
            N = 20;
            obj.kuka.goal_state = xG;
            
            % setup dircol
            t_min = 1;
            t_max = 6;
            options.integration_method = DirtranTrajectoryOptimization.MIDPOINT;
            prog = DirtranTrajectoryOptimization(obj.kuka, N, [t_min t_max],options);% arbitrary timescale
            
             % state constraint
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xG),N);
            
            % cost functions
            running_cost_func = @(dt,x,u) obj.kuka.runningCost(dt,x,u,xG);
            final_cost_func = @(dt,x) obj.kuka.finalCost(dt,x,xG);
            prog = prog.addRunningCost(running_cost_func);
            prog = prog.addFinalCost(final_cost_func);
            
            % joint limits
            [jl_min, jl_max] = obj.kuka.getJointLimits;
            x_min = [jl_min; -.75*ones(Nv,1)];
            x_max = [jl_max;  .75*ones(Nv,1)];
            prog = prog.addConstraint(BoundingBoxConstraint(repmat(x_min,N,1),repmat(x_max,N,1)), prog.x_inds);
            
            % additional solver options
            prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-3);
            prog.setSolverOptions('snopt','minorfeasibilitytolerance',1e-3);
            prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-3);
            
            % initial guess at the trajectory
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xG)]));

            % solve with dircol
            disp('Solving with dircol');
            tic;
            [xtraj, utraj, ~, ~, info] = prog.solveTraj(tf0,traj_init);
            toc
            disp('Finished dircol')
            
            % unpack the trajectory
            t = xtraj.getBreaks();
            x = xtraj.eval(xtraj.getBreaks());
            u = utraj.eval(utraj.getBreaks());
            
            % for debugging
            x
            info
            traj = xtraj.setOutputFrame(obj.kuka.getStateFrame());
            v = obj.kuka.constructVisualizer();
            playback(v,traj,struct('slider',true));

            
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