function runTrajOpt_vertical_throwing_non_robust_minimum_constraint

options=struct();
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = true;
options.ignore_self_collisions = true;
options.multiple_contacts = false;
options.active_collision_options.terrain_only = true;
 
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

global iteration_num
%global robustLCPcost_coeff

%% forward simulation
%trial 1, initial gripper pose is open
% q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.08; ...
%       0;0.79;0.09;0;0;0];
%trial 2, initial gripper pose is close
q0 = [-1.575;-1.4;0;1.27;0.0;1.1;0;0.06; ...
    0.0145;0.79;0.09;0;0;0];
%trial 5, inital gripper pose is open
% q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.06; ...
%       0.01;0.79;0.09;0;0;0];
x0 = [q0;zeros(nv,1)];
v.draw(0,x0);
kinematics_options.compute_gradients = 0;
kinsol = doKinematics(r, q0, [], kinematics_options);
iiwa_link_7_init = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
fr1 = r.forwardKin(kinsol,r.findLinkId('right_finger'),[0;0.04;0.1225],0);
R_ee = rpy2rotmat(iiwa_link_7_init(4:6));
rel_pos_object_gripper(1:3) = R_ee'*(q0(9:11) - iiwa_link_7_init(1:3));
rel_rot_object_gripper = rpy2rotmat(q0(12:14))*rpy2rotmat(iiwa_link_7_init(4:6));
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
q1 = q0;
q1(2) = q1(2) + 0.3;
%q1(1) = q0(1) + 0.8; 
q1(4) = q1(4) + 0.2; 
q1(6) = q1(6) - 0.1;
%q1(8) = q0(8) - 0.02;
kinsol = doKinematics(r, q1, [], kinematics_options);
iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
q1(9:11) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
q1(12:14) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
%trial 5
%q1 = q0;
%q1(8) = q1(8) - 0.015;
q1(11) = q1(11) - 0.03;% lift the height
x1 = [q1;zeros(nv,1)];
v.draw(0,x1);

qm_object = q1(9:14);
qm_object(3) = qm_object(3) + 0.08;

u0 = r.findTrim(q0);
u0(8) = 0;%-5;
u1 = r.findTrim(q1);
u1(8) = 0;%-5;
 
T0 = 2;
N = 15;
Nm = 7;

options.robustLCPcost_coeff = 1000;
options.Px_coeff = 0.1;
options.K = [10*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object),2*ones(nq_arm,nq_arm),zeros(nq_arm,nq_object)];
options.kappa = 1;
options.contact_robust_cost_coeff = 1e-8;
 
% ikoptions = IKoptions(r);
t_init = linspace(0,T0,N);
x_init = zeros(length(x0),N);
for i=1:length(x0)
    x_init(i,:) = linspace(x0(i,:),x1(i,:),N);
end
%run fwd IK for grasped object position
for i=2:N
    kinsol = doKinematics(r, x_init(:,i), [], kinematics_options);
    iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
    R_ee = rpy2rotmat(iiwa_link_7_final(4:6));
    x_init(9:11,i) = iiwa_link_7_final(1:3) + R_ee*rel_pos_object_gripper(1:3)';
    x_init(12:14,i) = rotmat2rpy((rel_rot_object_gripper*rpy2rotmat(iiwa_link_7_final(4:6))')');
end

traj_init.x = PPTrajectory(foh(t_init,x_init));
traj_init.x = traj_init.x.setOutputFrame(r.getStateFrame);

u_init = zeros(length(u0),N);
for i=1:length(u0)
    u_init(i,:) = linspace(u0(i,:),u1(i,:),N);
end
traj_init.u = PPTrajectory(foh(t_init,u_init));
traj_init.u = traj_init.u.setOutputFrame(r.getInputFrame);

%traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
%traj_init.u = PPTrajectory(zoh([0 T0],[u0, u0]));
T_span = [1 T0];

% Allcons = cell(0,1);
% [xtraj_init,snopt_info_ik,infeasible_constraint] = inverseKinTraj(r,t_init,traj_init.x,traj_init.x,ikoptions);
% v.playback(traj_init.x,struct('slider',true));

traj_opt = RobustContactImplicitTrajectoryOptimization_Kuka(r,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(q0_lb,q0_ub),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x1),N);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm_object),Nm,9:14);% constraint the object position and pose in the air during the middle phase
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(0.07),Nm,8);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(0.65,0.08),Nm,8);
%traj_opt = traj_opt.addStateConstraint(BoundingBoxConstraint(0.07,0.08),Nm+1,8);
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
%traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',100000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-4);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-4);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

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
global phi_cache_full 

    function [f,df] = running_cost_fun(h,x,u)
        R = 1e-6*eye(nu);
        Q = blkdiag(10*eye(8),0*eye(6),10*eye(14));
        g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
        f = h*g;
        df = [g, h*(x-x1)'*Q, h*u'*R];
    end

    function [f,df] = final_cost_fun(h,x)
        Qf = 1000*blkdiag(10*eye(8),0*eye(6),10*eye(14));
        g = (1/2)*(x-x1)'*Qf*(x-x1);
        f = h*g;
        df = [g, h*(x-x1)'*Qf];
    end

    function displayTraj(h,x,u,contact_force,LCP_slack_var)
        ts = [0;cumsum(h)];
        for i=1:length(ts)
            v.drawWrapper(0,x(:,i));
            pause(h(1)/3);
        end
        
        LCP_slack_var = LCP_slack_var';
        LCP_slack_var = [LCP_slack_var, LCP_slack_var(:,end)];
        fprintf('sum of slack variables along traj: %4.6f\n',sum(LCP_slack_var,2));
%         global robustLCPcost_coeff
%         if isempty(iteration_num)
%             robustLCPcost_coeff = 1;
%             iteration_num = 1;
%         elseif iteration_num > 15
%             robustLCPcost_coeff = 10;
%         elseif iteration_num > 30
%             robustLCPcost_coeff = 100;
%         elseif iteration_num > 45
%             robustLCPcost_coeff = 1000;    
%         end        
%         iteration_num = iteration_num + 1;
            
    end

end