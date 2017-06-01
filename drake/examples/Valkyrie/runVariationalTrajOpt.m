function [p,xtraj,utraj,z,F,info,traj_opt] = runVariationalTrajOpt()

robot_options = struct(); 

robot_options = applyDefaults(robot_options, struct('use_bullet', false,...
                                                    'terrain', RigidBodyFlatTerrain,...
                                                    'floating', true,...
                                                    'ignore_self_collisions', true,...
                                                    'ignore_friction', true,...
                                                    'enable_fastqp', false,...
                                                    'urdf_modifications_file', '',...
                                                    'dt', 0.001));
% silence some warnings
warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints')
warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits')

% construct robot model
p = Valkyrie(fullfile(getDrakePath,'examples','Valkyrie','urdf','urdf','valkyrie_simple.urdf'),robot_options);
p = p.removeCollisionGroupsExcept({'heel','toe'});
p = compile(p);

v = p.constructVisualizer;

nq = p.getNumPositions();
nv = p.getNumVelocities();
nx = nq+nv;
nu = p.getNumInputs();

% data=load(fullfile(getDrakePath,'examples','Valkyrie','data','valkyrie_fp_june2015_30joints_one_neck.mat'));

% Load nominal data
x0 = zeros(p.getNumStates,1);
x0(3) = 1.1734; 
x0(nq+1) = 0.1; 
q0 = x0(1:nq);
v0 = x0(nq+(1:nv));
T0 = 1.5;
N = 10;

% ----- Initial Guess ----- %
q1 = [0.4;q0(2:end)];
x1 = [q1;zeros(nv,1)];

t_init = linspace(0,T0,N);
% traj_init.x = PPTrajectory(foh([0 T0/2 T0],[x0, xm, x1]));
traj_init.x = PPTrajectory(foh([0 T0],[x0, x1]));
traj_init.u = PPTrajectory(zoh(t_init,0.1*randn(nu,N)));
T_span = [.25 T0];

options.add_ccost = true;
traj_opt = VariationalTrajectoryOptimization(p,N,T_span,options);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);  
% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qm),7);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q1),N);
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(v0),1);

[q_lb, q_ub] = getJointLimits(p);
q_ub(3) = q0(3) + 0.03;
q_lb(3) = q0(3) - 0.15;
traj_opt = traj_opt.addPositionConstraint(BoundingBoxConstraint(q_lb,q_ub),2:N-1);


state_cost = Point(getStateFrame(p),ones(nx,1));
state_cost.base_x = 0;
state_cost.base_y = 1;
state_cost.base_z = 5;
state_cost.base_pitch = 10;
state_cost.base_roll = 10;
state_cost.base_yaw = 10;
state_cost.leftAnklePitch=2;
state_cost.leftHipRoll = 2;
state_cost.leftHipPitch = 0;
state_cost.leftKneePitch = 0;
state_cost.rightAnklePitch=2;
state_cost.rightHipRoll = 2;
state_cost.rightHipPitch = 0;
state_cost.rightKneePitch = 0;
state_cost = double(state_cost);
Q = diag(state_cost); 


% traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(qf(1:6)),N);
%traj_opt = traj_opt.addInputConstraint(ConstantConstraint(zeros(3,1)),1:N-1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
% traj_opt = traj_opt.addPositionConstraint(periodic_constraint,{[1 N]});

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

save('val_traj.mat','xtraj','utraj','ctraj','btraj','straj','z','F');
v.playback(xtraj,struct('slider',true));

function [f,df] = running_cost_fun(h,x,u)
  R = 10*eye(nu);
  g = (1/2)*(x-x1)'*Q*(x-x1) + (1/2)*u'*R*u;
  f = h*g;
  df = [g, h*(x-x1)'*Q, h*u'*R];
end



  function displayTraj(h,x,u)
  
    ts = [0;cumsum(h)];
    for i=1:length(ts)
      v.drawWrapper(0,x(:,i));
      pause(h(1));
    end
   
end
end
