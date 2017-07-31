function RobustcontactImplicitBrick(visualize,position_tol,velocity_tol)
% tests that the contact implicit trajectory optimization can reproduce a
% simulation of the falling brick
rng(0)
if nargin < 1, visualize = false; end
if nargin < 2, position_tol = 1.5e-2; end
if nargin < 3, velocity_tol = 1e-1; end

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
plant = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
warning(w);
x0 = [0;0;1.0;0;0;0;1.5;zeros(5,1)];
%x0 = [0;0;1.0;0;0;0;zeros(6,1)];%free fall
xf = [1;0;0.5;0;0;0;zeros(6,1)];

N=30; tf=2.5;

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

scale_sequence = [10;1;.001;0];

for i=1:length(scale_sequence)
    scale = scale_sequence(i);
    
    options.compl_slack = scale*.01;
    options.lincompl_slack = scale*.001;
    options.jlcompl_slack = scale*.01;
    
    prog = RobustContactImplicitTrajectoryOptimization_Brick(plant,N,tf,options);
    %prog = ContactImplicitTrajectoryOptimization(plant,N,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','SuperbasicsLimit',1000);

    % prog = prog.setCheckGrad(true);
    
    %   snprint('snopt.out');
    
    % initial conditions constraint
    prog = addStateConstraint(prog,ConstantConstraint(x0),1);
    prog = addStateConstraint(prog,ConstantConstraint(xf),N);
    prog = prog.addTrajectoryDisplayFunction(@displayTraj);
    prog = prog.addRunningCost(@running_cost_fun);

    if i == 1,
        traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));
        traj_init.F_ext = PPTrajectory(foh([0,tf], 0.01*ones(1,2)));
        %traj_init.LCP_slack = PPTrajectory(foh([0,tf], 0.01*ones(1,2)));
    else
        traj_init.x = xtraj;
        traj_init.l = ltraj;
        traj_init.F_ext = F_exttraj;
    end
    [xtraj,utraj,ltraj,~,F_exttraj,z,F,info] = solveTraj(prog,tf,traj_init);
    
    if visualize
        v.playback(xtraj,struct('slider',true));
        % Create an animation movie
        %v.playbackAVI(xtraj, 'throwingBrick.avi');

        ts = getBreaks(xtraj);
        F_exttraj_data = F_exttraj.eval(ts);
        x_traj_data = xtraj.eval(ts);
        
    end
end

    function [f,df] = running_cost_fun(h,x,force)
        f = h*force'*force;
        df = [force'*force zeros(1,12) 2*h*force'];
    end

    function displayTraj(h,x,u,force)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(ts(i),x(:,i));
            xlim([-1.5, 6])
            pause(h(1)/5);
        end
        
%         LCP_slack = LCP_slack';
%         LCP_slack = [LCP_slack, LCP_slack(:,end)];
%         nominal_linewidth = 2.5;
%         color_line_type = 'r-';
%         figure(3)
%         plot(ts, LCP_slack(1,:), color_line_type, 'LineWidth',nominal_linewidth);
%         xlabel('t');
%         ylabel('slack variable');
%         hold off;
%         %         figure(4)
%         %         plot(ts, LCP_slack(2,:), color_line_type, 'LineWidth',nominal_linewidth);
%         %         xlabel('t');
%         %         ylabel('slack variable');
%         %         hold off;
%        fprintf('sum of slack variables along traj: %4.4f\n',sum(LCP_slack,2));
%        slack_sum_vec = [slack_sum_vec sum(LCP_slack,2)];
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