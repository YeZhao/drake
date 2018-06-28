function runTrajOpt2
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
% q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.08; ...
%       0;0.79;0.09;0;0;0];
%trial 2, initial gripper pose is close
% q0 = [-1.575;-1.4;0;1.27;0.0;1.1;0;0.057; ...
%     0.015;0.79;0.09;0;0;0];
%trial 3
qm = [-1.575;-1.4;0;1.27;0.0;1.1;0;0.06; ...
    0.0145;0.79;0.09;0;0;0];
%trial 5, inital gripper pose is open
% q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.06; ...
%       0.01;0.79;0.09;0;0;0];
xm = [qm;zeros(nv,1)];
v.draw(0,xm);
kinematics_options.compute_gradients = 0;
kinsol = doKinematics(r, qm, [], kinematics_options);
iiwa_link_7_init = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
fr1 = r.forwardKin(kinsol,r.findLinkId('right_finger'),[0;0.04;0.1225],0);
R_ee = rpy2rotmat(iiwa_link_7_init(4:6));
rel_pos_object_gripper(1:3) = R_ee'*(qm(9:11) - iiwa_link_7_init(1:3));
rel_rot_object_gripper = rpy2rotmat(qm(12:14))*rpy2rotmat(iiwa_link_7_init(4:6));
%xtraj_ts = simulate(r,[0 2],x0);
%v.playback(xtraj_ts,struct('slider',true));

%q1 = q0;
%q1(9) = q0(9)-0.1;
%q1(11) = q0(11)+0.1;
%q1(1:7) = q1(1:7) - [-1.57;-0.4;0;0.6;0;0.4;0];%[-1.57;-0.4;0;0.6;0;0.4;0];
%trial 1
% q1 = [0.7850;-0.6;0;1.27;0.0;0.35;0;0.08; ...
%       -0.57;-0.57;0.59;0;0;0];
%trial 2
%q1 = [-1.4;-1.4;0;1.27;0.0;1.1;0;0.06; ...
%      -0.124;0.78;0.09;0;0;0];
%trial 4
q1 = qm;
q1(2) = q1(2) + 0.3;
q1(1) = q1(1) + 0.8;
q1(6) = q1(6) - 0.25;
%q1(8) = q1(8) - 0.02;
kinsol = doKinematics(r, q1, [], kinematics_options);
iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
q1(9:11) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
q1(12:14) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
%trial 5
%q1 = q0;
%q1(8) = q1(8) - 0.015;
x1 = [q1;zeros(nv,1)];
%v.draw(0,x1);

%reposition initial state
q0 = qm;
q0(2) = q0(2) + 0.2;
q0(4) = q0(4) + 0.4;
q0(8) = 0.08;
x0 = [q0;zeros(nv,1)];
v.draw(0,x0);

u0 = r.findTrim(q0);
u0(8) = -5;
um = r.findTrim(qm);
um(8) = -5;
%u1 = r.findTrim(q1);
%u1(8) = -5;

T0 = 1;
N = 10;%phase 1: pick

r.uncertainty_source = '';%'+object_initial_position';%'object_initial_position'
if strcmp(r.uncertainty_source, 'friction_coeff') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_position')
    w_mu = load('friction_coeff_noise.dat');
    r.uncertain_mu_set = w_mu;
    r.uncertain_mu_mean = mean(r.uncertain_mu_set);
end
if strcmp(r.uncertainty_source, 'object_initial_position') || strcmp(r.uncertainty_source, 'friction_coeff+object_initial_position')
    w_phi = load('initial_position_noise.dat');
    r.uncertain_position_set = w_phi;
    r.uncertain_position_mean = mean(w_phi,2);
end
 
options.contact_robust_cost_coeff = 1;%important, if it is 0.1, can not solve successfully.
options.Px_coeff = 0.09;
options.Px_regularizer_coeff = 1e-1;
options.robustLCPcost_coeff = 1000;
options.K = [10*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object),2*sqrt(10)*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object)];
options.N = N;

% ikoptions = IKoptions(r);
t_init = linspace(0,T0,N);
x_init = zeros(length(x0),N);

%% phase 1
for i=1:length(x0)
    x_init(i,:) = linspace(x0(i,:),xm(i,:),N);
end

%run fwd kinematics for grasped object position
% for i=2:N1
%     kinsol = doKinematics(r, x_init1(:,i), [], kinematics_options);
%     iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
%     R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
%     x_init1(9:11,i) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
%     x_init1(12:14,i) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
% end
 
u_init = zeros(length(u0),N);
for i=1:length(u0)
    u_init(i,:) = linspace(u0(i,:),um(i,:),N);
end

traj_init.x = PPTrajectory(foh(t_init,x_init));
traj_init.x = traj_init.x.setOutputFrame(r.getStateFrame);
traj_init.u = PPTrajectory(foh(t_init,u_init));
traj_init.u = traj_init.u.setOutputFrame(r.getInputFrame);

% qm = q0;
% qm(2) = q0(2) + 0.4;
% qm(6) = q1(6) - 0.2;
% kinsol_m = doKinematics(r, qm, [], kinematics_options);
% iiwa_link_7_final_m = r.forwardKin(kinsol_m,r.findLinkId('iiwa_link_7'),[0;0;0],1);
% R_ee_m = rpy2rotmat(iiwa_link_7_final_m(4:6));
% qm(9:11) = iiwa_link_7_final_m(1:3) + R_ee_m*rel_pos_object_gripper(1:3)';
% qm(12:14) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final_m(4:6))')');
% qm(11) = qm(11) + 0.1;% adjust the height
% xm = [qm;zeros(nv,1)];
% v.draw(0,xm);

%traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
%traj_init.u = PPTrajectory(zoh([0 T0],[u0, u0]));
T_span = [1 T0];
% v.playback(traj_init.x,struct('slider',true));

% x0_ub = [q0;inf*ones(14,1)];
% x0_lb = [q0;-inf*ones(14,1)];
% x1_ub = [q1;inf*ones(14,1)];
% x1_lb = [q1;-inf*ones(14,1)];

% xfinal_lb = x1 - 0.05*ones(length(x1),1);
% xfinal_ub = x1 + 0.05*ones(length(x1),1);
xm_lb = xm - 0.1*ones(length(xm),1);
xm_ub = xm + 0.1*ones(length(xm),1);
xm_lb(8) = xm(8);
xm_ub(8) = xm(8);

traj_opt = RobustContactImplicitTrajectoryOptimization_Kuka(r,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(q0_lb,q0_ub),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x1),N);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(xm),N);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xm_lb,xm_ub),N);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xfinal_lb,xfinal_ub),N);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(xm_lb,xm_ub),N/2);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:7)),N,1:7);% free the finger final position
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(9:14)),N,9:14);
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
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',100000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',5e-3);%5e-5
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',5e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',5e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',5e-3);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);
global time_step
time_step = T0/(N-1); 

persistent sum_running_cost
persistent cost_index

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc
v.playback(xtraj,struct('slider',true));
keyboard

% % simulate with LQR gains
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

h_nominal = z(traj_opt.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z components
u_nominal = utraj.eval(t_nominal);
c_nominal = ctraj.eval(t_nominal);
c_normal_nominal = c_nominal(1:6:end,:);

    function [f,df] = running_cost_fun(h,x,u)
        R = 1e-6*eye(nu);
        Q = blkdiag(10*eye(8),0*eye(6),10*eye(14));
        g = (1/2)*(x-xm)'*Q*(x-xm) + (1/2)*u'*R*u;
        f = h*g;
        df = [g, h*(x-xm)'*Q, h*u'*R];
        
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
        Qf = 1000*blkdiag(10*eye(8),0*eye(6),10*eye(14));
        g = (1/2)*(x-xm)'*Qf*(x-xm);
        f = h*g;
        df = [g, h*(x-xm)'*Qf];
        fprintf('final cost: %4.4f\n',f);
        disp('------------------')
    end

    function displayTraj(h,x,u,contact_force,LCP_slack_var)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(h(1)/3); 
        end 
        
        if isempty(iteration_index)
            iteration_index = 1;
        else
            iteration_index = iteration_index + 1;
        end
        fprintf('iteration index: %4d\n',iteration_index);
        
        LCP_slack_var = LCP_slack_var';
        LCP_slack_var = [LCP_slack_var, LCP_slack_var(:,end)];
        if any(LCP_slack_var < 0)
            disp('here')
        end
          
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