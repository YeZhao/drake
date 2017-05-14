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
        
        function [utraj,xtraj]=swingUpTrajectory(obj)
            x0 = zeros(4,1);
            xf = double(obj.xG);
            tf0 = 4;%[changed]
            
            N = 21;
%            prog = DircolTrajectoryOptimization(obj,N,[2 6]);
           prog = DirtranTrajectoryOptimization(obj,N,[2 6]);
%             prog = RobustSamplingDirtranTrajectoryOptimization(obj,N,[2 6]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            prog = prog.addFinalCost(@finalCost);
            
            traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
            
            for attempts=1:10
                attempts
                tic
                size(traj_init)
                [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
                toc
                if info==1, break; end
            end
            
            % generate LQR gain matrix
            D = 2^2;
            E0 = zeros(4);
            Q = diag([10 10 1 1]);
            R = .1;
            Qf = 100*eye(4);
            obj = setInputLimits(obj,-10,10);
            prog_LQR = RobustSamplingDirtranTrajectoryOptimization(obj,N,D,E0,Q,R,Qf,[tf0 tf0]);
            
            %{reshape([obj.h_inds'; obj.x_inds(:,1:end-1); obj.u_inds],[],1); obj.x_inds(:,end)}
            
            h_vector = tf0/(N-1)*ones(1,N-1);
            t_span = linspace(0,tf0,N);
            xtraj_eval = xtraj.eval(t_span);
            utraj_eval = utraj.eval(t_span);
            state_full = reshape([h_vector;xtraj_eval(:,1:end-1);utraj_eval(:,1:end-1)],[],1);
            state_full = [state_full;xtraj_eval(:,end)];
            %zeros(nW,1)
            
            %state_full = 
            [K,~,~,~,E,dK,~,~,~,dE] = prog_LQR.deltaLQR(state_full,xf);
            
            function [g,dg] = cost(dt,x,u)
                persistent count;
                if isempty(count)
                    count = 1;
                else
                    count = count + 1;
                end
                count
                R = 1;
                g = sum((R*u).*u,1);
                dg = [zeros(1,1+size(x,1)),2*u'*R];
                return;
                
                xd = repmat([pi;0;0;0],1,size(x,2));
                xerr = x-xd;
                xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
                
                Q = diag([10,10,1,1]);
                R = 80;%100;
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
        
        function [utraj,xtraj,z,prog]=swingUpDirtran(obj,N)
            x0 = zeros(4,1); tf0 = 5; xf = double(obj.xG);
            
            obj = setInputLimits(obj,-10,10);
            prog = RobustDirtranTrajectoryOptimization(obj,N,1,zeros(4),eye(4),1,eye(4),[tf0 tf0]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            %prog = prog.addFinalCost(@finalCost);
            
            prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-4);
            prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-4);
            prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-4);
            
            function [g,dg] = cost(dt,x,u)
                Q = dt*1*eye(4);
                R = dt*.1;
                g = .5*(x-xf)'*Q*(x-xf) + .5*u'*R*u;
                dg = [g, (x-xf)'*Q, u'*R];
            end
            
            function [h,dh] = finalCost(t,x)
                Qf = 100*eye(4);
                h = .5*(x-xf)'*Qf*(x-xf);
                dh = [0, (x-xf)'*Qf];
            end
            
            traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
            tic
            [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
            toc
        end
        
        function [utraj,xtraj,z,prog]=robustSwingUp(obj,N,D,E0,Q,R,Qf)
            x0 = zeros(4,1); tf0 = 5; xf = double(obj.xG);
            
            obj = setInputLimits(obj,-10,10);
            prog = RobustDirtranTrajectoryOptimization(obj,N,D,E0,Q,R,Qf,[tf0 tf0]);
            prog = prog.addStateConstraint(ConstantConstraint(x0),1);
            prog = prog.addStateConstraint(ConstantConstraint(xf),N);
            prog = prog.addRunningCost(@cost);
            %prog = prog.addFinalCost(@finalCost);
            
            prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-3);
            prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-4);
            prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-4);
            
            prog = prog.addRobustCost(Q,1,Qf);
            prog = prog.addRobustInputConstraint();
            
            function [g,dg] = cost(dt,x,u)
                Q = dt*1*eye(4);
                R = dt*.1;
                g = .5*(x-xf)'*Q*(x-xf) + .5*u'*R*u;
                dg = [g, (x-xf)'*Q, u'*R];
            end
            
            function [h,dh] = finalCost(t,x)
                Qf = 100*eye(4);
                h = .5*(x-xf)'*Qf*(x-xf);
                dh = [0, (x-xf)'*Qf];
            end
            
            %traj_init.u = uguess;
            %traj_init.x = xguess;
            traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
            tic
            [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
            toc
        end
        
    end
    
end
