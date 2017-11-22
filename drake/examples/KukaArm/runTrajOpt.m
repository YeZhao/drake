function runTrajOpt

options=struct();
options.terrain = RigidBodyFlatTerrain();
options.use_bullet = true;
options.ignore_self_collisions = true;
options.multiple_contacts = false;
options.active_collision_options.terrain_only = false;

% options.with_weight = true;
% options.with_shelf_and_boxes = true;
r = KukaArm(options);

nq = r.getNumPositions();
nv = r.getNumVelocities();
nx = nq+nv;
nu = r.getNumInputs();

v=r.constructVisualizer;

%% forward simulation
%trial 1, initial gripper pose is open
% q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.08; ...
%       0;0.79;0.09;0;0;0];
%trial 2, initial gripper pose is close
q0 = [-1.57;-1.4;0;1.27;0.0;1.1;0;0.056; ...
      0.011;0.79;0.09;0;0;0];
x0 = [q0;zeros(nv,1)];
v.draw(0,x0);
kinematics_options.compute_gradients = 0;
kinsol = doKinematics(r, q0, [], kinematics_options);
iiwa_link_7_init = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
rel_pos_object_gripper = q0(9:14) - iiwa_link_7_init;
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
q1(2) = q0(2) + 0.02;
kinsol = doKinematics(r, q1, [], kinematics_options);
iiwa_link_7_final = r.forwardKin(kinsol,r.findLinkId('iiwa_link_7'),[0;0;0],1);
q1(9:14) = iiwa_link_7_final + rel_pos_object_gripper;
x1 = [q1;zeros(nv,1)];
%v.draw(0,x1);

u0 = r.findTrim(q0);
u0(8) = -10;

T0 = 2;
N = 8;

options.robustLCPcost_coeff = 100;

t_init = linspace(0,T0,N);
traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
traj_init.u = PPTrajectory(zoh([0 T0],[u0, u0]));
T_span = [1 T0];

traj_opt = RobustContactImplicitTrajectoryOptimization_Kuka(r,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addFinalCost(@final_cost_fun);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);  
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(x0),1);
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(1:7)),N,1:7);% free the finger final position
%traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1(9:14)),N,9:14);

[q_lb, q_ub] = getJointLimits(r);
% q_lb = max([q_lb, q0-0.2*ones(14,1)]')';
% q_ub = min([q_ub, q0+0.2*ones(14,1)]')';
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),1:N);
% u_ub = [inf;inf;inf;inf;inf;inf;inf;u0(8)];
% u_lb = [-inf;-inf;-inf;-inf;-inf;-inf;-inf;u0(8)];
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
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-3);
traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-3);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

v.playback(xtraj,struct('slider',true));

h_nominal = z(traj_opt.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z components
u_nominal = utraj.eval(t_nominal)';

function [f,df] = running_cost_fun(h,x,u)
  R = 1e-6*eye(nu);
  Q = blkdiag(1*eye(7),100,0*eye(6),10*eye(14));
  g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
  f = h*g;
  df = [g, h*(x-x1)'*Q, h*u'*R];
end
 
function [f,df] = final_cost_fun(h,x)
  Qf = 100*blkdiag(1*eye(8),0*eye(6),10*eye(14));
  g = (1/2)*(x-x1)'*Qf*(x-x1);
  f = h*g;
  df = [g, h*(x-x1)'*Qf];
end

function displayTraj(h,x,u,LCP_slack_var)
  ts = [0;cumsum(h)];
  for i=1:length(ts)
    v.drawWrapper(0,x(:,i));
    pause(h(1)/3);
  end
end

end