function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testNewTrajOpt_robustSampleTerrain(xtraj,utraj,ltraj,ljltraj)
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
options.terrain = RigidBodyStepTerrainVaryingHeight(0);
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);

paramstd = 1/5; % Standard deviation of the parameter value percent error
SampleNum = 1; % number of sampled terrain height
% perturb model parameters
paramerr = [];
for i = 1:SampleNum
    paramerr(i) = -0.05;%randn(1,1)*paramstd;
    if (paramerr(i) > 0.1)
        paramerr(i) = 0.1;
    elseif (paramerr(i) < -0.1)
        paramerr(i) = -0.1;
    end
    options(i).terrain = RigidBodyStepTerrainVaryingHeight(paramerr(i));
    options(i).floating = true;
    options(i).ignore_self_collisions = true;
    p_perturb(i) = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options(i));
    % trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);
end

v = p_perturb.constructVisualizer;
%v = p_perturb(1).constructVisualizer;
p = p_perturb;

%todo: add joint limits, periodicity constraint

N = 30;
T = 5;
T0 = 5;

% periodic constraint
R_periodic = zeros(p.getNumStates,2*p.getNumStates);
%R_periodic(2,2) = 1; %z
R_periodic(3,3) = 1; %pitch-hip w/symmetry
R_periodic(3,5) = 1; %pitch-hip w/symmetry
R_periodic(4,6) = 1; %knee w/symmetry
R_periodic(6,4) = 1; %knee w/symmetry
R_periodic(5,5) = -1; %hip w/symmetry

R_periodic(7,7) = 1; %x-vel
R_periodic(8,8) = 1; %z-vel
R_periodic(9,9) = 1; %pitch-hip w/symmetry
R_periodic(9,11) = 1; %pitch-hip w/symmetry
R_periodic(10,12) = 1; %knee w/symmetry
R_periodic(12,10) = 1; %knee w/symmetry
R_periodic(11,11) = -1; %hip w/symmetry

%R_periodic(2:end,p.getNumStates+2:end) = -eye(p.getNumStates-1);
R_periodic(3:end,p.getNumStates+3:end) = -eye(p.getNumStates-2);

periodic_constraint = LinearConstraint(zeros(p.getNumStates,1),zeros(p.getNumStates,1),R_periodic);

% x0 = [0;0;1;zeros(15,1)];
% xf = [0;0;1;zeros(15,1)];
x0 = [0;1;zeros(10,1)];
xf = [.4;1.;zeros(10,1)];

N2 = floor(N/2);

if nargin < 2
    %Try to come up with a reasonable trajectory
    %x1 = [.3;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
    x1 = [.3;1;-pi/8;pi/3;-pi/3;pi/3;zeros(6,1)];
    t_init = linspace(0,T0,N);
    %   traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
    traj_init.x = PPTrajectory(foh(t_init,[linspacevec(x0,x1,N2), linspacevec(x1,xf,N-N2)]));
    traj_init.u = PPTrajectory(foh(t_init,randn(3,N)));
    traj_init.l = PPTrajectory(foh(t_init,[repmat([1;zeros(7,1)],1,N2) repmat([zeros(4,1);1;zeros(3,1)],1,N-N2)]));
    traj_init.ljl = PPTrajectory(foh(t_init,zeros(p.getNumJointLimitConstraints,N)));
    traj_init.LCP_slack = PPTrajectory(foh(t_init, 0.01*ones(2,N)));
else
    t_init = xtraj.pp.breaks;
    traj_init.x = xtraj;
    traj_init.u = utraj;
    traj_init.l = ltraj;
    traj_init.ljl = ljltraj;
end
T_span = [1 T];

x0_min = [x0(1:5);-inf; -inf; 0; -inf(4,1)];
x0_max = [x0(1:5);inf;  inf; 0; inf(4,1)];
xf_min = [.4;-inf(11,1)];
xf_max = inf(12,1);
% xf_min = [xf(1:5);-inf(7,1)];
% xf_max = [xf(1:5);inf(7,1)];

scale = 0.1;
to_options.nlcc_mode = 2;
to_options.lincc_mode = 1;
to_options.compl_slack = scale*.01;
to_options.lincompl_slack = scale*.001;
to_options.jlcompl_slack = scale*.01;
to_options.lambda_mult = p.getMass*9.81*T0/N;
to_options.lambda_jl_mult = T0/N;

traj_opt = RobustContactImplicitTrajectoryOptimization(p,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xf),N);
traj_opt = traj_opt.addStateConstraint(periodic_constraint,{[1 N]});


traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);
slack_sum_vec = [];% vector storing the slack variable sum

% for i = 1:SampleNum
%     traj_opt = traj_opt.addRobustLCPConstraints(p_perturb(i),i,SampleNum);% [Ye: modify SampleNum parameter]
% end

% % extra robust cost function
% traj_opt.nonlincompl_constraints{i} = NonlinearComplementarityConstraint(@nonlincompl_fun,nX + traj_opt.nC,traj_opt.nC*(1+obj.nD),traj_opt.options.nlcc_mode,traj_opt.options.compl_slack);
% traj_opt.nonlincompl_slack_inds{i} = traj_opt.num_vars+1:traj_opt.num_vars + traj_opt.nonlincompl_constraints{i}.n_slack;
% traj_opt = traj_opt.addConstraint(traj_opt.nonlincompl_constraints{i},[traj_opt.x_inds(:,i+1);gamma_inds;lambda_inds]);

% traj_opt = traj_opt.setCheckGrad(true);
snprint('snopt.out');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
tic
[xtraj,utraj,ltraj,ljltraj,slacktraj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

v.playback(xtraj,struct('slider',true));
% Create an animation movie
%v.playbackAVI(xtraj, 'trial1.avi');

h_nominal = z(traj_opt.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z components
u_nominal = utraj.eval(t_nominal)';

% plot nominal model trajs
nominal_linewidth = 2.5;
color_line_type = 'b-';
figure(1)
subplot(3,1,1)
hold on;
plot(t_nominal', u_nominal(:,1), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('u1');
hold on;

subplot(3,1,2)
hold on;
plot(t_nominal', u_nominal(:,1), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('u2');
hold on;

subplot(3,1,3)
hold on;
plot(t_nominal', u_nominal(:,1), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('u3');
hold on;

figure(2)
subplot(2,3,1)
hold on;
plot(t_nominal', x_nominal(1,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('CoM x')
hold on;

subplot(2,3,2)
hold on;
plot(t_nominal', x_nominal(2,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('CoM z')
hold on;

subplot(2,3,3)
hold on;
plot(t_nominal', x_nominal(3,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Hip joint 1')
hold on;

subplot(2,3,4)
hold on;
plot(t_nominal', x_nominal(4,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 1')
hold on;

subplot(2,3,5)
hold on;
plot(t_nominal', x_nominal(5,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 1')
hold on;

subplot(2,3,6)
hold on;
plot(t_nominal', x_nominal(6,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 2')
hold on;

disp('finish traj opt')

    function [f,df] = running_cost_fun(h,x,u)
        f = h*u'*u;
        df = [u'*u zeros(1,12) 2*h*u'];
    end

    function displayTraj(h,x,u,LCP_slack)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(ts(i),x(:,i));
                  pause(h(1)/10);
        end
        
        LCP_slack = [LCP_slack, LCP_slack(:,end)];
        nominal_linewidth = 2.5;
        color_line_type = 'r-';
        figure(3)
        plot(ts, LCP_slack, color_line_type, 'LineWidth',nominal_linewidth);
        xlabel('t');
        ylabel('slack variable');
        hold off;
        fprintf('sum of slack variables along traj: %4.4f, %4.4f\n',sum(LCP_slack,2));
        slack_sum_vec = [slack_sum_vec sum(LCP_slack,2)];
    end
end