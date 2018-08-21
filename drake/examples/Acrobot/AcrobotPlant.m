classdef AcrobotPlant < Manipulator
    
    properties
        % parameters from Spong95 (except inertias are now relative to the
        % joints)
        % axis)
        l1 = 1; l2 = 2;
        m1 = 1; m2 = 1;
        g = 9.81;
        b1=.1;  b2=.1;
        %    b1=0; b2=0;
        lc1 = .5; lc2 = 1;
        Ic1 = .083;  Ic2 = .33;
        
        xG
        uG
        uncertainty_source
    end
    
    methods
        function obj = AcrobotPlant
            obj = obj@Manipulator(2,1);
            obj = setInputLimits(obj,-10,10);
            
            obj = setInputFrame(obj,CoordinateFrame('AcrobotInput',1,'u',{'tau'}));
            obj = setStateFrame(obj,CoordinateFrame('AcrobotState',4,'x',{'theta1','theta2','theta1dot','theta2dot'}));
            obj = setOutputFrame(obj,obj.getStateFrame);
            
            obj.xG = Point(obj.getStateFrame,[pi;0;0;0]);
            obj.uG = Point(obj.getInputFrame,0);
            
            %       obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',6,'p',...
            %         { 'b1','b2','lc1','lc2','Ic1','Ic2' }));
            %      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',4,'p',...
            %        {'m1','m2','Ic1','Ic2' }));
            obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',4,'p',...
                {'m1','m2','Ic1','Ic2'}));
            %    obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',1,'p',...
            %        {'m2' }));
            %      obj = setParamFrame(obj,CoordinateFrame('AcrobotParams',1,'p',...
            %        { 'lc2' }));
            obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
        end
        
        function [H,C,B] = manipulatorDynamics(obj,q,qd)
            % keep it readable:
            m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
            I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
            m2l1lc2 = m2*l1*lc2;  % occurs often!
            
            c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:));
            
            h12 = I2 + m2l1lc2*c(2);
            H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
            
            C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
            G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
            
            % accumulate total C and add a damping term:
            C = C*qd + G + [b1;b2].*qd;
            
            B = [0; 1];
        end
        
        function [T,U] = energy(obj,x)
            m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
            I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
            q = x(1:2); qd = x(3:4);
            c = cos(q(1:2,:));  s = sin(q(1:2,:));  c12 = cos(q(1,:)+q(2,:));
            
            T = .5*I1*qd(1)^2 + .5*(m2*l1^2 + I2 + 2*m2*l1*lc2*c(2))*qd(1)^2 + .5*I2*qd(2)^2 + (I2 + m2*l1*lc2*c(2))*qd(1)*qd(2);
            U = -m1*g*lc1*c(1) - m2*g*(l1*c(1)+lc2*c12);
        end
        
        % todo: also implement sodynamics here so that I can keep the
        % vectorized version?
        
        function [f,df,d2f,d3f] = dynamics(obj,t,x,u)
            f = dynamics@Manipulator(obj,t,x,u);
            if (nargout>1)
                [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
            end
        end
        
        function [f,df,d2f] = dynamics_w(obj,t,x,u,w)
            q = x(1:2);  qd = x(3:4);
            
            m1=obj.m1; m2=obj.m2; l1=obj.l1; g=obj.g; lc1=obj.lc1; lc2=obj.lc2; b1=obj.b1; b2=obj.b2;
            I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
            m2l1lc2 = m2*l1*lc2;  % occurs often!
            
            c = cos(q(1:2,:));  s = sin(q(1:2,:));  s12 = sin(q(1,:)+q(2,:));
            
            h12 = I2 + m2l1lc2*c(2);
            H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
            
            C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
            G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
            
            % accumulate total C and add a damping term:
            C = C*qd + G + [b1;b2].*qd;
            
            Bu = [0; 1];
            Bw = Bu;
            
            qdd = -H\(C - Bu*u - Bw*w);
            
            f = [qd; qdd];
            
            if (nargout>1)
                [df,d2f] = dynamicsGradients_w(obj,t,x,u,w,nargout-1);
            end
        end
        
        function nW = getNumDisturbances(obj)
            nW = 1;
        end
        
        function x = getInitialState(obj)
            x = zeros(4,1);%.1*randn(4,1);
        end
        
        function n = getNumPositions(obj)
            n = 2;
        end
        
        function n = getNumVelocities(obj)
            n = 2;
        end
        
        function [c,V]=balanceLQR(obj)
            Q = diag([10,10,1,1]); R = 1;
            if (nargout<2)
                c = tilqr(obj,obj.xG,obj.uG,Q,R);
            else
                if any(~isinf([obj.umin;obj.umax]))
                    error('currently, you must disable input limits to estimate the ROA');
                end
                [c,V] = tilqr(obj,obj.xG,obj.uG,Q,R);
                pp = feedback(obj.taylorApprox(0,obj.xG,obj.uG,3),c);
                options.method='levelSet';
                V=regionOfAttraction(pp,V,options);
            end
        end
        
        function [utraj,xtraj,z,prog]=swingUpTrajectory(obj,N)
            x0 = zeros(4,1);
            xf = double(obj.xG);
            tf0 = 6;
            
            prog = DirtranTrajectoryOptimization(obj,N,[6 6]);
            %prog = DircolTrajectoryOptimization(obj,N,[6 6]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);
            
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
            
            for attempts=1:1
                tic
                [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
                toc
                if info==1, break; end
            end
            
            function [g,dg] = cost(dt,x,u)
                R = 1;
                g = sum((R*u).*u,1);
                dg = [zeros(1,1+size(x,1)),2*u'*R];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Q = diag([10,10,1,1]);
                R = 100;
                g = sum((Q*xerr).*xerr + (R*u).*u,1);
                
                if (nargout>1)
                    dgddt = 0;
                    dgdx = 2*xerr'*Q;
                    dgdu = 2*u'*R;
                    dg = [dgddt,dgdx,dgdu];
                end
            end
            
            function [h,dh] = finalCost(t,x)
                h = t;
                dh = [1,zeros(1,size(x,1))];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Qf = 100*diag([10,10,1,1]);
                h = sum((Qf*xerr).*xerr,1);
                
                if (nargout>1)
                    dh = [0, 2*xerr'*Qf];
                end
            end
        end
        
        function [utraj,xtraj,z,prog,K]=robustSwingUpTrajectory_dirtran(obj,N,utrajArray,xtrajArray)
            x0 = zeros(4,1);
            xf = double(obj.xG);
            tf0 = 6;%[changed]
            
            nq = obj.getNumPositions;
            nv = obj.getNumVelocities;

            global sum_running_cost
            global sum_state_running_cost
            global sum_control_running_cost
            global cost_index
            global iteration_index
            iteration_index = 0;

            % LQR gains
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
             
            v = AcrobotVisualizer(obj);
             
            options.Px_coeff = 0.09;
            options.alpha = 0.5;
            options.kappa = 1;
            options.K = [0,10,0,sqrt(10)*2];
            options.contact_robust_cost_coeff = 1;%works with 0.5*randn noise.
            %tune alpha, contact_robust_cost_coeff, number of knot points.
             
            obj.uncertainty_source = 'physical_parameter_uncertainty';
            
            %prog = RobustDirtranTrajectoryOptimization(obj,N,[tf0 tf0],options);
            prog = RobustDircolTrajectoryOptimization(obj,N,Q,R,Qf,[tf0 tf0],options);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);
            prog = prog.addMotionDisplayFunction(@displayTraj);
            
            %prog = prog.setSolverOptions('snopt','MajorIterationsLimit',10000);
            %prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
            %prog = prog.setSolverOptions('snopt','IterationsLimit',100000000);
            %prog = prog.setSolverOptions('snopt','SuperbasicsLimit',1000000);
            %prog = prog.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
            %prog = prog.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);
            %prog = prog.setSolverOptions('snopt','MinorOptimalityTolerance',1e-3);
            %prog = prog.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);

            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
             
            warm_start = 0;
            if warm_start
                load('robust_test_alpha_p9_robust_cost_coeff_1_knot_point_30_warm_start.mat');
                traj_init.x = xtraj;%PPTrajectory(foh(t_init,x_nominal));
                traj_init.x = traj_init.x.setOutputFrame(obj.getStateFrame);
                
                traj_init.u = utraj;%PPTrajectory(foh(t_init,u_nominal));
                traj_init.u = traj_init.u.setOutputFrame(obj.getInputFrame);
                options.alpha = 0.8;
                v=AcrobotVisualizer(obj);
                iteration_index = 0;
                cost_index = [];
            end
 
            tic
            %size(traj_init)
            [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
            toc
            v.playback(xtraj,struct('slider',true));
            keyboard
            
            h = z(prog.h_inds);
            t = [0; cumsum(h)];
            x = xtraj.eval(t);
            u = utraj.eval(t)';
            
            % plot nominal model trajs
            nominal_linewidth = 2.5;
            color_line_type = 'b-';
            figure(10)
            hold on;
            plot(t', u, color_line_type, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('u');
            hold on;
            
            figure(22)
            subplot(2,2,1)
            hold on;
            plot(t', x(1,:), color_line_type, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('x_1')
            hold on;
            
            subplot(2,2,2)
            hold on;
            plot(t', x(2,:), color_line_type, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('x_2')
            hold on;
            
            subplot(2,2,3)
            hold on;
            plot(t', x(3,:), color_line_type, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('xdot_1')
            hold on;
            
            subplot(2,2,4)
            hold on;
            plot(t', x(4,:), color_line_type, 'LineWidth',nominal_linewidth);
            xlabel('t');
            ylabel('xdot_2')
            hold on;
            keyboard

            % %% pd control
            % Kp = options.K(2);
            % Kd = options.K(4);
            %
            % Kp = 500;
            % Kd = sqrt(Kp)*2;
            % sys=pdcontrol(obj,[Kp],[Kd],[2]);
            %
            % ltisys = pdcontrol(obj,Kp,Kd,[2]);
            % ltisys = setOutputFrame(ltisys,obj.getInputFrame);
            %
            % sys=feedback(obj,ltisys);
            
            %% LQR stabilization
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            
            %ytraj = simulate(cascade(setOutputFrame(xtraj, sys.getInputFrame), sys), [0, ts(end)], xstar_hand);
            
            num_trial = 100;
            num_success = 0;
            state_noise = load('random_state_noise.dat');
            param_uncertainty = load('physical_param_pertubation_simulate.dat');
            for i=1:num_trial
                obj.m1 = 1;
                obj.m2 = 1;
                obj.l1 = 1;
                obj.l2 = 2;
                obj.b1 = 0.1;
                obj.b2 = 0.1;
                obj.lc1 = 0.5;
                obj.lc2 = 1;
                obj.Ic1 = 0.0830;
                obj.Ic2 = 0.3300;
                 
                obj.m1 = obj.m1 + obj.m1*param_uncertainty(i,1)/10;
                obj.m2 = obj.m2 + obj.m2*param_uncertainty(i,2)/10;
                %obj.l1 = obj.l1 + obj.l1*param_uncertainty(i,3)/10;
                %obj.l2 = obj.l2 + obj.l2*param_uncertainty(i,4)/10;
                %obj.plant.b1  = obj.plant.b1 + obj.plant.b1*param_uncertainty(j,5);
                %obj.plant.b2  = obj.plant.b2 + obj.plant.b2*param_uncertainty(j,6);
                obj.lc1 = obj.lc1 + obj.lc1*param_uncertainty(i,7)/10;
                obj.lc2 = obj.lc2 + obj.lc2*param_uncertainty(i,8)/10;
                %obj.Ic1 = obj.Ic1 + obj.Ic1*param_uncertainty(i,9)/10;
                %obj.Ic2 = obj.Ic2 + obj.Ic2*param_uncertainty(i,10)/10;
                
                ltvsys = tvlqr(obj,xtraj,utraj,Q,R,Qf);
                sys=feedback(obj,ltvsys);
                
                xtraj_new = simulate(sys,xtraj.tspan, [0;0;0;0]+state_noise(:,i));%
                
                %xtraj_new = simulate(cascade(setOutputFrame(xtraj, sys.getInputFrame), sys),xtraj.tspan, [0;0;0;0]+state_noise(:,i));
                
                v.playback(xtraj_new,struct('slider',true));
                ts = getBreaks(xtraj_new);
                xtraj_new_data = xtraj_new.eval(ts);
                if sum(abs(xtraj_new_data(:,end) - xf)) < 1e-2
                    num_success = num_success + 1;
                end
            end
            keyboard
            
            % N = 31, direct collocation
            %if alpha = 1, num_success = 57; mean deviation 72.1, covariance = 43.17
            %if alpha = 0.9 and no warm start, num_success = 59; mean deviation = 54.1, covariance = 51.7
            %if alpha = 0.8 and no warm start, num_success = 60; mean deviation = 62, covariance = 83.9, show 32 major interation limit
            %if alpha = 0.7 and no warm start, num_success = 66; mean deviation = 46, covariance = 155, not working yet
            
            % N = 61, direct collocation
            %if alpha = 1 and no warm start, num_success = ; mean deviation = 34, covariance = 10
            %if alpha = 0.9 and no warm start, num_success = 50; mean deviation = 35, covariance = 11.3
            %if alpha = 0.8 and no warm start, num_success = 49; mean deviation = 38, covariance = 12.5
            %if alpha = 0.7 and no warm start, num_success = 59; mean deviation = 35, covariance = 15
            
            %if using direct transcription, N = 31, alpha = 1, num_success = 4;
            %if using direct transcription, N = 101, alpha = 1, num_success = 29;
            
            %if warm start, alpha = 0.8, num_success = 30;
            
            %keyboard
            
            %% manually created PD stabilization
            % ts = getBreaks(xtraj);
            % h = tf0/(N-1);
            % utraj_data = utraj.eval(ts);
            % xtraj_data = xtraj.eval(ts);
            %
            % kp = 10;
            % kd = sqrt(kp)*2;
            %
            % K = [0,kp,0,kd];
            % g = 9.81;
            %
            % q_real(:,1) = xtraj_data(1:nq,1);
            % qdot_real(:,1) = xtraj_data(nq+1:nq+nv,1);
            % x_real_full(:,1) = xtraj_data(:,1);
            %
            % stabilitation_scenario = 'initial_state_noise';
            %
            % if isempty(stabilitation_scenario)
            %     sample_length = 1;
            % end
            %
            % if strcmp(stabilitation_scenario, 'initial_state_noise')
            %     w_x_init = load('initial_state_noise.dat');
            %     sample_length = size(w_x_init,2);
            % end
            %
            % for m=1:sample_length
            %     if strcmp(stabilitation_scenario, 'initial_state_noise')
            %         xtraj_data(:,1) = xtraj_data(:,1) + w_x_init(:,m);
            %         q_real(:,1) = xtraj_data(1:nq,1);
            %         qdot_real(:,1) = xtraj_data(nq+1:nq+nv,1);
            %     end
            %
            %     for i=1:N-1
            %         %feedback ctrl in q position
            %         F_fb(:,i) = K(:,1:nq)*(xtraj_data(1:nq,i) - q_real(:,i)) + K(:,nq+1:nq+nv)*(xtraj_data(1+nq:nq+nv,i) - qdot_real(:,i));
            %         %feedforward
            %         F_ff(:,i) = utraj_data(:,i);
            %         F_net(:,i) = F_fb(:,i) + F_ff(:,i);
            %
            %         [xdot,~] = obj.dynamics(0,x_real_full(:,i),F_net(:,i));
            %         xdn = x_real_full(:,i) + tf0/(N-1)*xdot;
            %
            %         x_real_full(:,i+1) = xdn;
            %         q_real(1:nq,i+1) = x_real_full(1:nq,i+1);
            %         qdot_real(1:nv,i+1) = x_real_full(1+nq:nq+nv,i+1);
            %     end
            %
            %     x_simulated = x_real_full;
            %     xtraj_simulated = PPTrajectory(foh(ts,x_simulated));
            %     xtraj_simulated = xtraj_simulated.setOutputFrame(obj.getStateFrame);
            %     v.playback(xtraj_simulated,struct('slider',true));
            %
            %     disp('-------------')
            %     fprintf('sample index: %4d\n',m);
            %     fprintf('final state deviation: %4.8f\n',norm(x_real_full(:,end) - xtraj_data(:,end)));
            %     fprintf('full trajectory state deviation cost: %4.8f\n',norm(x_real_full - xtraj_data));
            %     fprintf('final object height deviation cost: %4.8f\n',x_real_full(11,end) - xtraj_data(11,end));
            % end
            %keyboard
            
            % generate LQR gain matrix             
            % manually refactor the decision variable vector for deltaLQR input
            h_vector = tf0/(N-1)*ones(1,N-1);
            t_span = linspace(0,tf0,N);
            xtraj_eval = xtraj.eval(t_span);
            utraj_eval = utraj.eval(t_span);
            state_full = reshape([h_vector;xtraj_eval(:,1:end-1);utraj_eval(:,1:end-1)],[],1);
            %state_full = [state_full;xtraj_eval(:,end)];
            
            D = 2^2;
            E0 = zeros(4);
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            
            %if (nargout>4)
            %    % solve LQR feedback gain matrix of nominal model
            %    K = prog.deltaLQR(state_full,xf);
            %end
            
            function [g,dg] = cost(dt,x,u)
                % persistent count;
                % if isempty(count)
                %     count = 1;
                % else
                %     count = count + 1;
                % end
                % count
                
                %penalize control only
                R = 1;
                g = sum((R*u).*u,1);
                dg = [zeros(1,1+size(x,1)),2*u'*R];
                
                if isempty(cost_index)
                    cost_index = 1;
                    sum_running_cost = g;
                    %sum_state_running_cost = sum((Q*xerr).*xerr,1);
                    %sum_control_running_cost = sum((R*u).*u,1);
                elseif cost_index == N-2
                    sum_running_cost = sum_running_cost + g;
                    %sum_state_running_cost = sum_state_running_cost + sum((Q*xerr).*xerr,1);
                    %sum_control_running_cost = sum_control_running_cost + sum((R*u).*u,1);
                    fprintf('sum of running cost: %4.4f\n',sum_running_cost);
                    fprintf('sum of state running cost: %4.4f\n',sum_state_running_cost);
                    fprintf('sum of control running cost: %4.4f\n',sum_control_running_cost);
                    cost_index = [];
                else
                    sum_running_cost = sum_running_cost + g;
                    %sum_state_running_cost = sum_state_running_cost + sum((Q*xerr).*xerr,1);
                    %sum_control_running_cost = sum_control_running_cost + sum((R*u).*u,1);
                    cost_index = cost_index + 1;
                end
                
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Q = diag([1,1,0.1,0.1]);
                R = 0.1;
                g = sum((Q*xerr).*xerr,1) + sum((R*u).*u,1);
                
                if g ~= sum((Q*xerr).*xerr,1) + sum((R*u).*u,1);
                    keyboard
                end
                
                if (nargout>1)
                    dgddt = 0;
                    dgdx = 2*xerr'*Q;
                    dgdu = 2*u'*R;
                    dg = [dgddt,dgdx,dgdu];
                end
                
                %g_numeric = g;
                %dg_numeric = dg;
                %X0 = [dt,x',u];
                
                %[g_numeric,dg_numeric] = geval(@(X0) cost_check(X0),X0,struct('grad_method','numerical'));
                % valuecheck(dg,dg_numeric,1e-5);
                function [g,dg] = cost_check(X0)
                    dt = X0(1);
                    x = X0(2:5)';
                    u = X0(6);
                    
                    xd = repmat([pi;0;0;0],1,size(x,2));
                    xerr = x-xd;
                    xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                    
                    Q = diag([1,1,0.1,0.1]);
                    R = 0.1;
                    g = sum((Q*xerr).*xerr) + sum((R*u).*u,1);
                                       
                    if (nargout>1)
                        dgddt = 0;
                        dgdx = 2*xerr'*Q;
                        dgdu = 2*u'*R;
                        dg = [dgddt,dgdx,dgdu];
                    end
                    
                end
            end
            
            function [h,dh] = finalCost(t,x)
                h = t;
                dh = [1,zeros(1,size(x,1))];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Qf = 100*diag([10,10,1,1]);
                h = sum((Qf*xerr).*xerr,1);
                
                if (nargout>1)
                    dh = [0, 2*xerr'*Qf];
                end
            end
            
            function displayTraj(h,x,u)
                ts = [0;cumsum(h)];
                for i=1:length(ts)
                    v.drawWrapper(0,x(:,i));
                    %pause(h(1)/3);
                end
                
                iteration_index = iteration_index + 1;
                disp('------------------')
                fprintf('iteration index: %4d\n',iteration_index);
                
                % global robustLCPcost_coeff
                % if isempty(iteration_num)
                %     robustLCPcost_coeff = 1;
                %     iteration_num = 1;
                % elseif iteration_num > 15
                %     robustLCPcost_coeff = 10;
                % elseif iteration_num > 30
                %     robustLCPcost_coeff = 100;
                % elseif iteration_num > 45
                %     robustLCPcost_coeff = 1000;
                % end
                % iteration_num = iteration_num + 1;
            end
        end
        
        function [utraj,xtraj,z,prog,K]=swingUpTrajectory_dircol(obj,N,utrajArray,xtrajArray)
            x0 = zeros(4,1);
            xf = double(obj.xG);
            tf0 = 6;%[changed]
            
            % LQR gains
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            
            prog = RobustDircolTrajectoryOptimization(obj,N,Q,R,Qf,[tf0 tf0]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);
            
            if nargin > 2 % this part is used for re-generating nominal trajs
                Qr = diag([10 10 1 1]);
                Rr = .1;
                Qrf = 100*eye(4);
                prog = prog.addRobustAverageRunningCost(utrajArray,xtrajArray,Qr,Qrf,Rr);
            end
            
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
            
            for attempts=1:1
                attempts
                tic
                %size(traj_init)
                [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
                toc
                if info==1, break; end
            end
            
            % generate LQR gain matrix             
            % manually refactor the decision variable vector for deltaLQR input
            h_vector = tf0/(N-1)*ones(1,N-1);
            t_span = linspace(0,tf0,N);
            xtraj_eval = xtraj.eval(t_span);
            utraj_eval = utraj.eval(t_span);
            state_full = reshape([h_vector;xtraj_eval(:,1:end-1);utraj_eval(:,1:end-1)],[],1);
            %state_full = [state_full;xtraj_eval(:,end)];
            
            D = 2^2;
            E0 = zeros(4);
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            
            if (nargout>4)
                % solve LQR feedback gain matrix of nominal model
                K = prog.deltaLQR(state_full,xf);
            end
            
            function [g,dg] = cost(dt,x,u)
%                 persistent count;
%                 if isempty(count)
%                     count = 1;
%                 else
%                     count = count + 1;
%                 end
%                 count
                R = 1;
                g = sum((R*u).*u,1);
                dg = [zeros(1,1+size(x,1)),2*u'*R];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Q = diag([10,10,1,1]);
                R = 100;
                g = sum((Q*xerr).*xerr + (R*u).*u,1);
                
                if (nargout>1)
                    dgddt = 0;
                    dgdx = 2*xerr'*Q;
                    dgdu = 2*u'*R;
                    dg = [dgddt,dgdx,dgdu];
                end
            end
            
            function [h,dh] = finalCost(t,x)
                h = t;
                dh = [1,zeros(1,size(x,1))];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Qf = 100*diag([10,10,1,1]);
                h = sum((Qf*xerr).*xerr,1);
                
                if (nargout>1)
                    dh = [0, 2*xerr'*Qf];
                end
            end            
        end
        
        function [utraj,xtraj,z,prog]=robustSwingUpTrajectory_dircol(obj,utraj,xtraj,K,Qr,Qrf,Rr,N)
            x0 = zeros(4,1); tf0 = 6; xf = double(obj.xG);
                        
            % LQR gains
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            
            %obj = setInputLimits(obj,-10,10);
            prog = RobustDircolTrajectoryOptimization(obj,N,Q,R,Qf,[tf0 tf0]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);%[Ye: why this is commented out]
            
            % double check the terminal condition
            % prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-3);
            % prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-4);
            % prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-4);
            
            prog = prog.addRobustRunningCost(utraj,xtraj,K,Qr,Qrf,Rr);
            %prog = prog.addRobustInputConstraint();
            
            traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
            tic
            [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
            toc
            
            function [g,dg] = cost(dt,x,u)
                R = 1;
                g = sum((R*u).*u,1);
                dg = [zeros(1,1+size(x,1)),2*u'*R];
            end
            
            function [h,dh] = finalCost(t,x)
                h = t;
                dh = [1,zeros(1,size(x,1))];
            end
        end
    end
end
