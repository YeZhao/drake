function RobustcontactImplicitBrick(visualize,position_tol,velocity_tol)
% tests that the contact implicit trajectory optimization can reproduce a
% simulation of the falling brick
% rng(0)
if nargin < 1, visualize = false; end
if nargin < 2, position_tol = 1.5e-2; end
if nargin < 3, velocity_tol = 1e-1; end

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);
x0 = [0;0;2.0;0;0;0;.5;zeros(5,1)];
%x0 = [0;0;1.0;0;0;0;zeros(6,1)];%free fall
xf = [1;0;0.5;0;0;0;zeros(6,1)];

N=5; tf=2;

plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
xtraj_ts = simulate(plant_ts,[0 tf],x0);
x0 = xtraj_ts.eval(0);
warning(w);
if visualize
    v = constructVisualizer(plant_ts);
    v.playback(xtraj_ts);
end

options = struct();
options.integration_method = RobustContactImplicitTrajectoryOptimization_Brick.MIXED;
%options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

% scale_sequence = [1;.1;.01;.001;0];

% for i=1:length(scale_sequence)
%     scale = scale_sequence(i);
%     
%     options.compl_slack = scale*.01;
%     options.lincompl_slack = scale*.001;
%     options.jlcompl_slack = scale*.01;
    
    prog = RobustContactImplicitTrajectoryOptimization_Brick(plant_ts,N,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',20000);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',2000000);
    prog = prog.setSolverOptions('snopt','SuperbasicsLimit',10000);
    prog = prog.setSolverOptions('snopt','MajorOptimalityTolerance',1e-6);
    prog = prog.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-6);
    prog = prog.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-6);

    % prog = prog.setCheckGrad(true);
    
    %   snprint('snopt.out');
    
    % initial conditions constraint
    prog = addStateConstraint(prog,ConstantConstraint(x0),1);
    prog = addStateConstraint(prog,ConstantConstraint(xf),N);
    prog = prog.addTrajectoryDisplayFunction(@displayTraj);
    prog = prog.addRunningCost(@running_cost_fun);
    
    %     if i == 1,
    traj_init.x = PPTrajectory(foh([0,tf],[x0,xf]));
    traj_init.F_ext = PPTrajectory(foh([0,tf], 0.01*ones(1,2)));
    traj_init.LCP_slack = PPTrajectory(foh([0,tf], 0.01*ones(1,2)));
    slack_sum_vec = [];% vector storing the slack variable sum

    %     else
    %         traj_init.x = xtraj;
    %         traj_init.l = ltraj;
    %         traj_init.F_ext = F_exttraj;
    %     end
    [xtraj,utraj,ltraj,slacktraj,F_exttraj,z,F,info,infeasible_constraint_name] = solveTraj(prog,tf,traj_init);
    
    if visualize
        v.playback(xtraj,struct('slider',true));
        % Create an animation movie
        %v.playbackAVI(xtraj, 'throwingBrick.avi');

        ts = getBreaks(xtraj);
        F_exttraj_data = F_exttraj.eval(ts);
        xtraj_data = xtraj.eval(ts);
        ltraj_data = ltraj.eval(ts);
        nD = 4;
        nC = 8;
        lambda_n_data = ltraj_data(1:nD+2:end,:);
        
        %ltraj_data.
        figure(1)
        plot(ts, xtraj_data(3,:),'b-');
        hold on;
        for i=1:nC
            plot(ts, lambda_n_data(i,:),'r-');
            hold on;
        end
        
        figure(2)
        plot(ts, F_exttraj_data,'r-');
    end
% end

    function [f,df] = running_cost_fun(h,x,force)
        f = h*force'*force;
        df = [force'*force zeros(1,12) 2*h*force'];
        
        f_numeric = f;
        df_numeric = df;
%        disp('check gradient')
%         [f_numeric,df_numeric] = geval(@(h,x,force) running_cost_fun_check(h,x,force),h,x,force,struct('grad_method','numerical'));
%         valuecheck(df,df_numeric,1e-3);
%         valuecheck(f,f_numeric,1e-3);
        
        function [f,df] = running_cost_fun_check(h,x,force)
            f = h*force'*force;
            df = [force'*force zeros(1,12) 2*h*force'];
        end
    end

    function displayTraj(h,x,u,force,LCP_slack)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(ts(i),x(:,i));
            xlim([-1.5, 6])
            pause(h(1)/5);
        end
        
        LCP_slack = LCP_slack';
        LCP_slack = [LCP_slack, LCP_slack(:,end)];
        nominal_linewidth = 2.5;
        color_line_type = 'r-';
%         figure(3)
%         plot(ts, LCP_slack(1,:), color_line_type, 'LineWidth',nominal_linewidth);
%         xlabel('t');
%         ylabel('slack variable');
%         hold off;
%         
%         figure(4)
%         plot(ts, force, color_line_type, 'LineWidth',nominal_linewidth);
%         xlabel('t');
%         ylabel('external force');
%         hold off;
        
       fprintf('sum of slack variables along traj: %4.4f\n',sum(LCP_slack,2));
       fprintf('sum of external force along traj: %4.4f\n',sum(force));

       slack_sum_vec = [slack_sum_vec sum(LCP_slack,2)];
    end

% check if the two simulations did the same thing:
ts = getBreaks(xtraj_ts);
valuecheck(ts,getBreaks(xtraj));
xtraj_data = xtraj.eval(ts);
xtraj_ts_data = xtraj_ts.eval(ts);
nq = plant.getNumPositions();
nv = plant.getNumVelocities();
valuecheck(xtraj_data(1:nq,:),xtraj_ts_data(1:nq,:),position_tol); % is there a correct tolerance here?
valuecheck(xtraj_data(nq+(1:nv),:),xtraj_ts_data(nq+(1:nv),:),velocity_tol); % is there a correct tolerance here?

% TIMEOUT 750
end