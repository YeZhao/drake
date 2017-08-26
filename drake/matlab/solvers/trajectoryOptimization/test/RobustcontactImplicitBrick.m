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

N=30; tf=2;

%% instantiate RigidBodyTerrain with different heights
w_phi = load('terrain_height_noise5.dat');
%w_phi = normrnd(zeros(1,n_sig_point),sqrt(Pw(1,1)),1,n_sig_point);%height noise
%save -ascii terrain_height_noise5.dat w_phi
n_sig_point = 28;
for i=1:n_sig_point
    sample_options.terrain = RigidBodyFlatTerrain(w_phi(i));
    sample_options.floating = true;
    w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
    plant_sample{i} = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),sample_options);
    warning(w);
    plant.plant_sample{i} = plant_sample{i};% add multiple RigidBodyManipulators with Sampled Terrain Height into the normal RigidBodyManipulator
end

%% previous setting(August-22-17)
% x0 = [0;0;2.0;0;0;0;0.5;zeros(5,1)];
% %x0 = [0;0;1.0;0;0;0;zeros(6,1)];%free fall
% xf = [1;0;0.5;0;0;0;zeros(6,1)];%previous setting(August-22-17)
% xf_min = [1;0;-inf;0;0;0;zeros(6,1)];
% xf_max = [1;0;inf;0;0;0;zeros(6,1)];

%new setting(August-22-17)
x0 = [0;0;2.0;0;0;0;10;zeros(5,1)];
xf = [6.256;0;0.5;0;0;0;zeros(6,1)];
xf_min = [6.256;0;0.5;0;0;0;zeros(6,1)];
xf_max = [6.256;0;0.5;0;0;0;zeros(6,1)];

plant_ts = TimeSteppingRigidBodyManipulator(plant,tf/(N-1));
% w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
% xtraj_ts = simulate(plant_ts,[0 tf],x0);
% x0 = xtraj_ts.eval(0);
% warning(w);
if visualize
    v = constructVisualizer(plant_ts);
    %v.playback(xtraj_ts,struct('slider',true));
    
%     ts = getBreaks(xtraj_ts);
%     xtraj_ts_data = xtraj_ts.eval(ts);
%     
%     %ltraj_data.
%     figure(1)
%     plot(ts, xtraj_ts_data(3,:),'b--');
%     hold on;
%     plot(ts, xtraj_ts_data(1,:),'k--');
%     xlabel('t [s]','fontsize',15);ylabel('position/force','fontsize',15);
%     
%     figure(2)
%     plot(xtraj_ts_data(1,:), xtraj_ts_data(3,:),'b--');
%     xlabel('x [m]');ylabel('z [m]');
end

options = struct();
options.integration_method = RobustContactImplicitTrajectoryOptimization_Brick.MIXED;
%options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

options.contact_robust_cost_coeff = 0.0001;
options.robustLCPcost_coeff = 1000;
options.Px_coeff = 10;
options.Kx_gain = 5;
options.Kxd_gain = 5;
options.Kz_gain = 5;
options.Kzd_gain = 5;
options.kappa = 1;

persistent sum_running_cost
persistent cost_index

prog = RobustContactImplicitTrajectoryOptimization_Brick(plant_ts,N,tf,options);
prog = prog.setSolverOptions('snopt','MajorIterationsLimit',20000);
prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
prog = prog.setSolverOptions('snopt','IterationsLimit',2000000);
prog = prog.setSolverOptions('snopt','SuperbasicsLimit',10000);
prog = prog.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
prog = prog.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
prog = prog.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
%prog = prog.setCheckGrad(true);
 
%snprint('snopt.out');

% initial conditions constraint
prog = addStateConstraint(prog,ConstantConstraint(x0),1);
prog = addStateConstraint(prog,BoundingBoxConstraint(xf_min,xf_max),N);
prog = prog.addTrajectoryDisplayFunction(@displayTraj);
prog = prog.addRunningCost(@running_cost_fun);

traj_init.x = PPTrajectory(foh([0,tf],[x0,xf]));
traj_init.F_ext = PPTrajectory(foh([0,tf], 0.01*ones(2,2)));
traj_init.LCP_slack = PPTrajectory(foh([0,tf], 0.01*ones(1,2)));
slack_sum_vec = [];% vector storing the slack variable sum
 
[xtraj,utraj,ltraj,~,slacktraj,F_exttraj,z,F,info,infeasible_constraint_name] = solveTraj(prog,tf,traj_init);
 
if visualize
    v.playback(xtraj,struct('slider',true));
    % Create an animation movie
    %v.playbackAVI(xtraj, 'throwingBrick.avi');
    
    ts = getBreaks(xtraj);
    F_exttraj_data = F_exttraj.eval(ts)*h;%convert impulse to force.
    xtraj_data = xtraj.eval(ts);
    ltraj_data = ltraj.eval(ts);
    nD = 4;
    nC = 8;
    lambda_n_data = ltraj_data(1:nD+2:end,:);
    
    %ltraj_data.
    figure(1)
    hold on;
    plot(ts, xtraj_data(3,:),'b-');
    hold on;
    plot(ts, xtraj_data(1,:),'k-');
    hold on;
    for i=1:nC
        plot(ts, lambda_n_data(i,:),'r-');
        hold on;
    end
    legend('passive case z position','passive case x position','robust case z position','robust case x position','robust case ground reaction force');
    title('Time profile of positions and forces','fontsize',22);
    
    figure(2)
    hold on;
    plot(xtraj_data(1,:), xtraj_data(3,:),'b-');
    xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
    title('2D Cartesian CoM trajectory','fontsize',22)
    legend('passive case','robust case')
    ylim([0,2.1])
    
    figure(3)
    plot(ts, F_exttraj_data(1,:),'b-');
    hold on;
    plot(ts, F_exttraj_data(2,:),'r-');
    xlabel('t [s]','fontsize',20);ylabel('force [N]','fontsize',20);
    title('Force Profile','fontsize',22)
    legend('horizontal force','vertical force')
    ylim([-10.5,10.5])
end
    function [f,df] = running_cost_fun(h,x,force)
        cost_coeff = 1;
        
        f = cost_coeff*h*force'*force;
        df = [cost_coeff*force'*force zeros(1,12) 2*cost_coeff*h*force'];
        
        f_numeric = f;
        df_numeric = df;
        % disp('check gradient')
        % [f_numeric,df_numeric] = geval(@(h,x,force) running_cost_fun_check(h,x,force),h,x,force,struct('grad_method','numerical'));
        % valuecheck(df,df_numeric,1e-3);
        % valuecheck(f,f_numeric,1e-3);
        
        if isempty(cost_index)
            cost_index = 1;
            sum_running_cost = f;
        elseif cost_index == N-2
            sum_running_cost = sum_running_cost + f;
            fprintf('sum of running cost: %4.4f\n',sum_running_cost);
            disp('------------------')
            cost_index = [];
        else
            sum_running_cost = sum_running_cost + f;
            cost_index = cost_index + 1;
        end
        
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
        % figure(3)
        % plot(ts, LCP_slack(1,:), color_line_type, 'LineWidth',nominal_linewidth);
        % xlabel('t');
        % ylabel('slack variable');
        % hold off;
        %
        % figure(4)
        % plot(ts, force, color_line_type, 'LineWidth',nominal_linewidth);
        % xlabel('t');
        % ylabel('external force');
        % hold off;
        
        fprintf('sum of slack variables along traj: %4.4f\n',sum(LCP_slack,2));
        fprintf('sum of external x force along traj: %4.4f\n',sum(abs(force(1,:))));
        fprintf('sum of external z force along traj: %4.4f\n',sum(abs(force(2,:))));
        
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