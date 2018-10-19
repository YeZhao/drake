function runTrajOpt_objflip_internalDyn
options=struct();
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = true;
options.ignore_self_collisions = true;
options.multiple_contacts = false;
options.active_collision_options.terrain_only = true;

global iteration_index
iteration_index = 0;
global example_name;
example_name = 'kuka_arm';

% options.with_weight = true;
% options.with_shelf_and_boxes = true;
r = KukaArm(options);

nq = r.getNumPositions();
nv = r.getNumVelocities();
nx = nq+nv;
nu = r.getNumInputs();
nq_arm = 8;
nq_object = nq - nq_arm;

v=r.constructVisualizer;

%% forward simulation
%trial 1, initial gripper pose is open
q0 = [-1.575;-.93;0;1.57;0.0;-0.62;0;0.06; ...
    0.0145;0.58;0.06;0;0;0;0];
%q0 = [zeros(14,1);0];
x0 = [q0;zeros(nv,1)];
v.draw(0,x0);

kinematics_options.compute_gradients = 0;
kinsol = doKinematics(r, q0, [], kinematics_options);
iiwa_link_7_init = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
fr1 = r.forwardKin(kinsol,r.findLinkId('right_finger'),[0;0.04;0.1225],0);
R_ee = rpy2rotmat(iiwa_link_7_init(4:6));
rel_pos_object_gripper(1:3) = R_ee'*(q0(9:11) - iiwa_link_7_init(1:3));
rel_rot_object_gripper = rpy2rotmat(q0(12:14))*rpy2rotmat(iiwa_link_7_init(4:6));

qm = q0;
qm(6) = qm(6) + 1.82;
%q1(8) = q1(8) - 0.02;
kinsol = doKinematics(r, qm, [], kinematics_options);
iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
qm(9:11) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
qm(12:14) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
xm = [qm;zeros(nv,1)];
v.draw(0,xm);

q1 = [-1.575;-.93;0;1.57;0.0;1.6;0;0.06; ...
    0.0145;1.1;0.06;3.14159;0;0;0];
x1 = [q1;zeros(nv,1)];
v.draw(0,x1);

u0 = r.findTrim(q0);
u0(8) = -5;
um = r.findTrim(qm);
um(8) = -5;
u1 = r.findTrim(q1);
u1(8) = -5;

T0 = 2;
N = 35;%10;
N1 = 28;%phase 1: pick
N2 = N - N1;%phase 2: throw

r.uncertainty_source = '';%'friction_coeff+object_initial_position';%'object_initial_position'
r.uncertainty_source_default = r.uncertainty_source;
if strcmp(r.uncertainty_source, 'friction_coeff') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_position') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_orientation')
    w_mu = load('friction_coeff_noise.dat');
    r.uncertain_mu_set = w_mu;
    r.uncertain_mu_mean = mean(r.uncertain_mu_set);
end
if strcmp(r.uncertainty_source, 'object_initial_position') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_position')
    w_phi = load('initial_position_noise.dat');
    phi_scaling = 2;
    r.uncertain_position_set = w_phi/phi_scaling;
    r.uncertain_position_mean = mean(w_phi/phi_scaling,2);
end
if strcmp(r.uncertainty_source, 'object_initial_orientation') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_orientation')
    w_ori = load('initial_orientation_noise.dat');
    r.uncertain_orientation_set = w_ori;
    r.uncertain_orientation_mean = mean(w_ori,2);
end

options.contact_robust_cost_coeff = 0.1;%important, if it is 0.1, can not solve successfully.
options.Px_coeff = 0.09;
options.Px_regularizer_coeff = 1e-1;
options.robustLCPcost_coeff = 1000;
options.K = [10*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object),2*sqrt(10)*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object)];
options.N1 = N1;
options.test_name = 'pick_and_throw_motion';
options.alpha = 0.2;
options.kappa = 1;
options.robust_cost_type = 'ML cost';%1:non-ML cost; %2: ML cost

% ikoptions = IKoptions(r);
t_init = linspace(0,T0,N);
x_init = zeros(length(x0),N);

%% phase 1
for i=1:length(x0)
    x_init1(i,:) = linspace(x0(i,:),xm(i,:),N1);
end
%run fwd kinematics for grasped object position
for i=2:N1
    kinsol = doKinematics(r, x_init1(:,i), [], kinematics_options);
    iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
    R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
    x_init1(9:11,i) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
    x_init1(12:14,i) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
end
 
u_init1 = zeros(length(u0),N1);
for i=1:length(u0)
    u_init1(i,:) = linspace(u0(i,:),um(i,:),N1);
end
 
%% phase 2
for i=1:length(xm)
    x_init2(i,:) = linspace(xm(i,:),x1(i,:),N2);
end

u_init2 = zeros(length(um),N2);
for i=1:length(um)
    u_init2(i,:) = linspace(um(i,:),u1(i,:),N2);
end
 
x_init = [x_init1,x_init2];
u_init = [u_init1,u_init2];
traj_init.x = PPTrajectory(foh(t_init,x_init));
traj_init.x = traj_init.x.setOutputFrame(r.getStateFrame);
traj_init.u = PPTrajectory(foh(t_init,u_init));
traj_init.u = traj_init.u.setOutputFrame(r.getInputFrame);
 
%traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
%traj_init.u = PPTrajectory(zoh([0 T0],[u0, u0]));
T_span = [T0 10];
% v.playback(traj_init.x,struct('slider',true));

warm_start = 0;
if warm_start
    load('robust_test23_pick_and_place_motion_phi_scaling_1_good_motion_quadratic_robust_cost.mat');
    traj_init.x = PPTrajectory(foh(t_init,x_nominal));
    traj_init.x = traj_init.x.setOutputFrame(r.getStateFrame);

    traj_init.u = PPTrajectory(foh(t_init,u_nominal));
    traj_init.u = traj_init.u.setOutputFrame(r.getInputFrame);
    options.alpha = 1;
    options.test_name = 'pick_and_place_motion';
    v=r.constructVisualizer;
    iteration_index = [];
end

% x0_ub = [q0;inf*ones(14,1)];
% x0_lb = [q0;-inf*ones(14,1)];
% x1_ub = [q1;inf*ones(14,1)];
% x1_lb = [q1;-inf*ones(14,1)];

% xfinal_lb = x1 - 0.05*ones(length(x1),1);
% xfinal_ub = x1 + 0.05*ones(length(x1),1);
% xm_lb = xm - 0.05*ones(length(xm),1);
% xm_ub = xm + 0.05*ones(length(xm),1);
 
traj_opt = RobustContactImplicitTrajectoryOptimization_Kuka(r,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
%traj_opt = traj_opt.addFinalCost(@final_cost_fun2);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(q0_lb,q0_ub),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x1),N);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xm),N1);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(x0-0.01*ones(length(x0),1),x0+0.01*ones(length(x0),1)),1);
traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xm-0.01*ones(length(xm),1),xm+0.01*ones(length(xm),1)),N1);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xfinal_lb,xfinal_ub),N);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xm_lb,xm_ub),N/2);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:7)),N,1:7);% free the finger final position
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(9:14)),N,9:14);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(8:14)),N,8:14);

[q_lb, q_ub] = getJointLimits(r);
% q_lb = max([q_lb, q0-0.2*ones(14,1)]')';
% q_ub = min([q_ub, q0+0.2*ones(14,1)]')';
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),1:N);
% u_ub = [200*ones(7,1);-10];
% u_lb = [-200*ones(8,1)];
% traj_opt = traj_opt.addInputConstraint(BoundingBoxConstraint(u_lb,u_ub),1:N-1);

% ub_N = q1;
% ub_N(1:8) = q_ub(1:8);
% lb_N = q1;
% lb_N(1:8) = q_lb(1:8);
% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(lb_N,ub_N),N);
% traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(zeros(nv,1)),1);
% traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),2:N-1);

% state_cost = Point(getStateFrame(r),ones(nx,1));
% state_cost.base_x = 10;
% state_cost.base_y = 10;
% state_cost.base_z = 10;
% state_cost.base_pitch = 1;
% state_cost.base_roll = 1;
% state_cost.base_yaw = 1;
% state_cost = double(state_cost);
% Q = diag(state_cost); 

% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodi c_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',100000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','ScaleOption',0);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

if ~warm_start
    persistent sum_running_cost
    persistent cost_index
end

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc
v.playback(xtraj,struct('slider',true));
h_nominal = z(traj_opt.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z components
u_nominal = utraj.eval(t_nominal);
c_nominal = ctraj.eval(t_nominal);
c_normal_nominal = c_nominal(1:6:end,:);
lambda_n_data = [];
figure(1)
plot(t_nominal,x_nominal(8,:));
title('gripper position')
keyboard

%% stabilization
ts = getBreaks(xtraj);
h = T0/(N-1);
utraj_data = utraj.eval(ts);
xtraj_data = xtraj.eval(ts);
ctraj_data = ctraj.eval(ts);
nD = 4;
nC = 12;

kp = 500;
kd = sqrt(kp)*2;
kp_gripper = 1000;
kd_gripper = sqrt(kp_gripper)*2;

K = [blkdiag(kp*eye(nq_arm-1),kp_gripper),zeros(nq_arm,nq_object),blkdiag(kd*eye(nq_arm-1),kd_gripper),zeros(nq_arm,nq_object)];

g = 9.81;

q_real(:,1) = xtraj_data(1:nq,1);
qdot_real(:,1) = xtraj_data(nq+1:nq+nv,1);
x_real_full(:,1) = xtraj_data(:,1);

stabilitation_scenario = 'friction_coeff';

if isempty(stabilitation_scenario)
    sample_length = 1;
end

if strcmp(stabilitation_scenario, 'friction_coeff')
    w_mu = load('friction_coeff_noise.dat');
    sample_length = length(w_mu);
end

if strcmp(stabilitation_scenario, 'object_initial_position')
    w_phi = load('initial_position_noise.dat');
    %w_phi = w_phi/terrain_height_scale_factor;
    r.uncertain_position_set = w_phi;
    r.uncertain_position_mean = mean(w_phi,2);
    sample_length = length(w_phi);
end

if strcmp(stabilitation_scenario, 'friction_coeff+object_initial_position')
    w_mu = load('friction_coeff_noise.dat');
    w_phi = load('initial_position_noise.dat');
    %w_phi = w_phi/terrain_height_scale_factor;
    r.uncertain_position_set = w_phi;
    r.uncertain_position_mean = mean(w_phi);
    sample_length = length(w_phi);
end

for m=1:sample_length
    if strcmp(stabilitation_scenario, 'friction_coeff')
        r.uncertainty_source = 'friction_coeff';
        r.uncertain_mu = w_mu(m);
    elseif strcmp(stabilitation_scenario, 'object_initial_position')
        r.uncertainty_source = 'object_initial_position';
        r.uncertain_phi = w_phi(m,:);
    elseif strcmp(stabilitation_scenario, 'friction_coeff+object_initial_position')
        r.uncertainty_source = 'friction_coeff+object_initial_position';
        r.uncertain_phi = w_phi(m,:);
    end
    
    for i=1:N-1
        %feedback ctrl in q position
        F_fb(:,i) = K(:,1:nq)*(xtraj_data(1:nq,i) - q_real(:,i)) + K(:,nq+1:nq+nv)*(xtraj_data(1+nq:nq+nv,i) - qdot_real(:,i));
        %feedforward
        F_ff(:,i) = utraj_data(:,i);
        F_net(:,i) = F_fb(:,i) + F_ff(:,i);
        
        [xdn,~] = r.update(T0/(N-1),x_real_full(:,i),F_net(:,i));
        x_real_full(:,i+1) = xdn;
        q_real(1:nq,i+1) = x_real_full(1:nq,i+1);
        qdot_real(1:nv,i+1) = x_real_full(1+nq:nq+nv,i+1);
    end
     
    x_simulated = x_real_full;
    xtraj_simulated = PPTrajectory(foh(ts,x_simulated));
    xtraj_simulated = xtraj_simulated.setOutputFrame(r.getStateFrame);
    v.playback(xtraj_simulated,struct('slider',true));
    
    % figure(4)
    % hold on;
    % plot(x_real_full(1,:), x_real_full(3,:),'r-');
    % hold on;
    % plot(xtraj_data(1,:), xtraj_data(3,:),'b-');
    % xlabel('x [m]','fontsize',20);ylabel('z [m]','fontsize',20);
    % title('joint trajectory','fontsize',22)
    % %legend('passive case','robust case')
    % ylim([0,2.1])
    
    disp('-------------')
    fprintf('sample index: %4d\n',m);
    fprintf('final state deviation: %4.8f\n',norm(x_real_full(:,end) - xtraj_data(:,end)));
    fprintf('full trajectory state deviation cost: %4.8f\n',norm(x_real_full - xtraj_data));
    fprintf('final object height deviation cost: %4.8f\n',x_real_full(11,end) - xtraj_data(11,end));
end
keyboard

df = [];

%% simulate with LQR gains
% % LQR Cost Matrices
Q = diag(10*ones(1,nx));
R = .1*eye(nu);
Qf = 100*eye(nx);

ltvsys = tvlqr(r,xtraj,utraj,Q,R,Qf);
sys=feedback(r,ltvsys);
xtraj_new = simulate(sys,xtraj.tspan, x0);
v.playback(xtraj_new,struct('slider',true));
 
%% pd-control LTI trial
kp = 100;
kd = sqrt(kp)*1.5;

K = [kp*eye(nq_arm),kp*eye(nq_arm,nq_object),kd*eye(nq_arm),kd*eye(nq_arm,nq_object)];

ltisys = LinearSystem([],[],[],[],[],-K);
ltisys = setInputFrame(ltisys,CoordinateFrame([r.getStateFrame.name,' - ', mat2str(x0,3)],length(x0),r.getStateFrame.prefix));
r.getStateFrame.addTransform(AffineTransform(r.getStateFrame,ltisys.getInputFrame,eye(length(x0)),-x0));
ltisys.getInputFrame.addTransform(AffineTransform(ltisys.getInputFrame,r.getStateFrame,eye(length(x0)),+x0));
ltisys = setOutputFrame(ltisys,r.getInputFrame);

sys = feedback(r,ltisys);
% Forward simulate dynamics with visulazation, then playback at realtime
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
%sys = mimoCascade(sys,v,[],[],output_select);
warning(S);
xtraj_new = simulate(sys,xtraj.tspan,x0);
playback(v,xtraj_new,struct('slider',true));

    function [f,df] = running_cost_fun(h,x,u)
        R = 1e-6*eye(nu);
        Q = blkdiag(10*eye(8),0*eye(7),10*eye(15));
        g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
        f = h*g;
        df = [g, h*(x-x1)'*Q, h*u'*R];
        
        if isempty(cost_index)
            cost_index = 1;
            sum_running_cost = f;
        elseif cost_index == N-2
            sum_running_cost = sum_running_cost + f;
            fprintf('sum of running cost: %4.4f\n',sum_running_cost);
            cost_index = [];
        else
            sum_running_cost = sum_running_cost + f;
            cost_index = cost_index + 1;
        end
    end

    function [f,df] = final_cost_fun(h,x)
        Qf = 1000*blkdiag(10*eye(8),0*eye(7),10*eye(15));
        g = (1/2)*(x-x1)'*Qf*(x-x1);
        f = h*g;
        df = [g, h*(x-x1)'*Qf];
        fprintf('final cost: %4.4f\n',f);
        disp('------------------')
    end

    function [f,df] = final_cost_fun2(h,x)
        x_far = 10*x1;
        Qf = 1000*blkdiag(10*eye(8),0*eye(6),10*eye(14));
        g = (1/2)*(x-x_far)'*Qf*(x-x_far);
        f = h*g;
        df = [g, h*(x-x_far)'*Qf];
        fprintf('final cost: %4.4f\n',f);
        disp('------------------')
    end

    function displayTraj(h,x,u,contact_force,LCP_slack_var)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(h(1)/3);
        end
        
        LCP_slack_var = LCP_slack_var';
        LCP_slack_var = [LCP_slack_var, LCP_slack_var(:,end)];
        
        iteration_index = iteration_index + 1;
        fprintf('iteration index: %4d\n',iteration_index);
         
        fprintf('sum of slack variables along traj: %4.6f\n',sum(LCP_slack_var,2));
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