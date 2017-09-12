function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = testNewTrajOpt_robustSampleTerrain_new_multiple_steps(xtraj,utraj,ltraj,ljltraj)
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
    paramerr(i) = 0.15;%randn(1,1)*paramstd;
    if (paramerr(i) > 0.3)
        paramerr(i) = 0.3;
    elseif (paramerr(i) < -0.3)
        paramerr(i) = -0.3;
    end
    options(i).terrain = RigidBodyStepTerrainMultipleSteps();%RigidBodyStepTerrainVaryingHeight(paramerr(i));
    options(i).floating = true;
    options(i).ignore_self_collisions = true;
    p_perturb(i) = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options(i));
    % trajopt = ContactImplicitTrajectoryOptimization(p,[],[],[],10,[1 1]);
end 

v = p_perturb.constructVisualizer;
%v = p_perturb(1).constructVisualizer;
p = p_perturb;
nq = p.getNumPositions(); 
nv = p.getNumVelocities();
nx = nq+nv;
nu = p.getNumInputs();

w_phi = load('terrain_height_noise5.dat');
%w_phi = normrnd(zeros(1,n_sig_point),sqrt(Pw(1,1)),1,n_sig_point);%height noise
%save -ascii terrain_height_noise5.dat w_phi 
n_sig_point = 28;
for i=1:n_sig_point
    sample_options.terrain = RigidBodyStepTerrainVaryingHeight(w_phi(i));
    sample_options.floating = true;
    w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
    p_sample{i} = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',sample_options);
    warning(w);
    p.plant_sample{i} = p_sample{i};% add multiple RigidBodyManipulators with Sampled Terrain Height into the normal RigidBodyManipulator
end

%todo: add joint limits, periodicity constraint

N = 70;
T = 3;
T0 = 3;
 
% periodic constraint
R_periodic = zeros(p.getNumStates,2*p.getNumStates);
R_periodic(2,2) = 1; %z
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

R_periodic(2:end,p.getNumStates+2:end) = -eye(p.getNumStates-1);
%R_periodic(3:end,p.getNumStates+3:end) = -eye(p.getNumStates-2);

periodic_constraint = LinearConstraint(zeros(p.getNumStates,1),zeros(p.getNumStates,1),R_periodic);

x0 = [0;1;zeros(10,1)];
xf = [1.2;1;zeros(10,1)];
%xf = [3.2;1.;zeros(10,1)];

N2 = floor(N/2);
 
if nargin < 2
    %Try to come up with a reasonable trajectory
    x1 = [0.6;1;pi/8-pi/16;pi/8;-pi/8;pi/8;zeros(6,1)];
    %x1 = [2;1;pi/8;pi/5;-pi/5;pi/5;zeros(6,1)];
    t_init = linspace(0,T0,N);
    %   traj_init.x = PPTrajectory(foh(t_init,linspacevec(x0,xf,N)));
    traj_init.x = PPTrajectory(foh(t_init,[linspacevec(x0,x1,N2), linspacevec(x1,xf,N-N2)]));
    traj_init.u = PPTrajectory(foh(t_init,zeros(3,N)));%randn(3,N)
    traj_init.l = PPTrajectory(foh(t_init,[repmat([1;zeros(7,1)],1,N2) repmat([zeros(4,1);1;zeros(3,1)],1,N-N2)]));
    traj_init.ljl = PPTrajectory(foh(t_init,zeros(p.getNumJointLimitConstraints,N)));
    traj_init.LCP_slack = PPTrajectory(foh(t_init, 0.01*ones(1,N)));
else
    t_init = xtraj.pp.breaks;
    traj_init.x = xtraj;
    traj_init.u = utraj;
    traj_init.l = ltraj;
    traj_init.ljl = ljltraj;
end
T_span = [1 T];

x0_min = [x0(1:5);-inf; 0; 0; -inf(4,1)];
x0_max = [x0(1:5);inf;  0; 0; inf(4,1)];
xf_min = [1;-inf(11,1)];
%xf_min = [3.2;-inf(11,1)];
xf_max = inf(12,1);

p_ts = TimeSteppingRigidBodyManipulator(p,T/(N-1));
w = warning('off','Drake:TimeSteppingRigidBodyManipulator:ResolvingLCP');
  
scale = 0.01;
to_options.nlcc_mode = 2;% robust mode %original: 2;
to_options.lincc_mode = 1; 
to_options.compl_slack = scale*.01;
to_options.lincompl_slack = scale*.001;
to_options.jlcompl_slack = scale*.01;
to_options.lambda_mult = p.getMass*9.81*T0/N;
to_options.lambda_jl_mult = T0/N;
 
to_options.contact_robust_cost_coeff = 1e-30;%1e-5;%0.0001;  
to_options.robustLCPcost_coeff = 1000;
to_options.Px_coeff = 0.01; 
%to_options.K = [zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,3)];%[three rows: hip,knee,knee]
to_options.K = [zeros(3,3),0.1*ones(3,3),zeros(3,3),0.1*ones(3,3)];%[three rows: hip,knee,knee]
to_options.kappa = 1;
running_cost_coeff = 1; 
to_options.add_ccost = true;
 
persistent sum_running_cost
persistent running_cost_index
persistent sum_foot_height_cost
persistent foot_height_cost_index

traj_opt = RobustContactImplicitTrajectoryOptimization(p_ts,N,T_span,to_options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addRunningCost(@foot_height_fun); 
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0_min,x0_max),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xf_min,xf_max),N);
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
 
%traj_opt = traj_opt.setCheckGrad(true); 
snprint('snopt.out'); 
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',500);%20000
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',20000);%200000
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',500000);%1000000
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',10000); 
traj_opt = traj_opt.setSolverOptions('snopt','VerifyLevel',0);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);

tic 
[xtraj,utraj,ltraj,ljltraj,slacktraj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc
snprint('snopt.out');

v.playback(xtraj,struct('slider',true));
xlim([-1.5, 6])
% Create an animation movie
%v.playbackAVI(xtraj, 'trial6_small_terrain_perturb_with_two_ERM_cost_full_nonlcompl_three_walking_constraints.avi');

h_nominal = z(traj_opt.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z components
u_nominal = utraj.eval(t_nominal)';
f_nominal = ltraj.eval(t_nominal);
slack_nominal = slacktraj.eval(t_nominal)';

%convert impulse into force
h = t_nominal(2)-t_nominal(1);
f_nominal = f_nominal./h;

fprintf('sum of x value: %4.4f\n',sum(sum(abs(x_nominal))));

foot_height = compute_foot_height(x_nominal);
fprintf('sum of left foot height integration value: %4.4f\n',sum(sum(abs(foot_height(1,:)))));
fprintf('sum of right foot height integration value: %4.4f\n',sum(sum(abs(foot_height(2,:)))));

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
ylabel('Hip joint 2')
hold on;

subplot(2,3,6)
hold on;
plot(t_nominal', x_nominal(6,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 2')
hold on;

figure(3)
subplot(2,3,1)
hold on;
plot(t_nominal', x_nominal(7,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('CoM x')
hold on;

subplot(2,3,2)
hold on;
plot(t_nominal', x_nominal(8,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('CoM z')
hold on;

subplot(2,3,3)
hold on;
plot(t_nominal', x_nominal(9,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Hip joint 1')
hold on;

subplot(2,3,4)
hold on;
plot(t_nominal', x_nominal(10,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 1')
hold on;

subplot(2,3,5)
hold on;
plot(t_nominal', x_nominal(11,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Hip joint 2')
hold on;

subplot(2,3,6)
hold on;
plot(t_nominal', x_nominal(12,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('Knee joint 2')
hold on;

figure(4)
subplot(2,2,1)
hold on;
plot(t_nominal', f_nominal(1,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f normal (left leg)')
hold on;

subplot(2,2,2)
hold on;
plot(t_nominal', f_nominal(2,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f tangential1 (left leg)')
hold on;

subplot(2,2,3)
hold on;
plot(t_nominal', f_nominal(3,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f tangential 2 (left leg)')
hold on;

subplot(2,2,4)
hold on;
plot(t_nominal', f_nominal(4,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('tangential velocity (left leg)')
hold on;

figure(5)
subplot(2,2,1)
hold on;
plot(t_nominal', f_nominal(5,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f normal (right leg)')
hold on;

subplot(2,2,2)
hold on;
plot(t_nominal', f_nominal(6,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f tangential1 (right leg)')
hold on;

subplot(2,2,3)
hold on;
plot(t_nominal', f_nominal(7,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('f tangential 2 (right leg)')
hold on;

subplot(2,2,4)
hold on;
plot(t_nominal', f_nominal(8,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('tangential velocity (right leg)')
hold on;
 
figure(6)
subplot(2,1,1)
hold on;
plot(t_nominal', foot_height(1,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('left foot height')
hold on;

subplot(2,1,2)
hold on;
plot(t_nominal', foot_height(2,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('right foot height')
hold on;

% simulate with LQR gains
% LQR Cost Matrices
Q = diag(10*ones(1,12));
R = .1*eye(3);
Qf = 100*eye(12); 

ltvsys = tvlqr(p,xtraj,utraj,Q,R,Qf);

sys=feedback(p,ltvsys);

xtraj_new = simulate(sys,xtraj.tspan, x0);
v.playback(xtraj_new,struct('slider',true));

disp('finish traj opt')

    function foot_height = compute_foot_height(x)
        q_full = x(1:nq,:);
        
        for i=1:size(q_full,2)
            [phi,~,~,~,~,~,~,~,n] = p.contactConstraints(q_full(:,i),false,struct('terrain_only',true));
            foot_height(:,i) = phi;
        end
        
    end
     
    function [f,df] = foot_height_fun(h,x,u)
        q = x(1:nq);
         
        [phi,~,~,~,~,~,~,~,n] = p.contactConstraints(q,false,struct('terrain_only',true));
        phi0 = [.2;.2];
        K = 30; 
        I = find(phi < phi0);
        f = K*(phi(I) - phi0(I))'*(phi(I) - phi0(I));
        % phi: 2x1
        % n: 2xnq
        df = [0 2*K*(phi(I)-phi0(I))'*n(I,:) zeros(1,nv+nu)];
        
        if isempty(foot_height_cost_index)
            foot_height_cost_index = 1;
            sum_foot_height_cost = f;
        elseif foot_height_cost_index == N-2
            sum_foot_height_cost = sum_foot_height_cost + f;
            fprintf('sum of foot height cost: %4.4f\n',sum_foot_height_cost);
            disp('------------------')
            foot_height_cost_index = [];
        else
            sum_foot_height_cost = sum_foot_height_cost + f;
            foot_height_cost_index = foot_height_cost_index + 1;
        end
        
        %    K = 100;
        %    K_log = 100;
        %    f = sum(-K*log(K_log*phi + .2));
        %    df = [0 sum(-K*K_log*n./(K_log*repmat(phi,1,length(q)) + .2)) zeros(1,15)];
    end

    function [f,df] = running_cost_fun(h,x,u)
        f = h*u'*u;
        df = [u'*u zeros(1,12) 2*h*u'];
        
        if isempty(running_cost_index)
            running_cost_index = 1;
            sum_running_cost = f;
        elseif running_cost_index == N-2
            sum_running_cost = sum_running_cost + f;
            fprintf('sum of running cost: %4.4f\n',sum_running_cost);
            running_cost_index = [];
        else
            sum_running_cost = sum_running_cost + f;
            running_cost_index = running_cost_index + 1;
        end
    end

    function displayTraj(h,x,u,LCP_slack)
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
        figure(1)
        plot(ts,x(1,:));
        xlabel('t');ylabel('x');
        
        figure(3)
        plot(ts, LCP_slack(1,:), color_line_type, 'LineWidth',nominal_linewidth);
        xlabel('t');
        ylabel('slack variable');
        hold off;
        %         figure(4)
        %         plot(ts, LCP_slack(2,:), color_line_type, 'LineWidth',nominal_linewidth);
        %         xlabel('t');
        %         ylabel('slack variable');
        %         hold off;
        fprintf('sum of slack variables along traj: %4.4f\n',sum(LCP_slack,2));
        slack_sum_vec = [slack_sum_vec sum(LCP_slack,2)];
    end
end