function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = variationalTrajOpt()

[~,xtraj,utraj] = testNewTrajOpt();

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);

N = 15;
T0 = 3;

% ----- Periodicity Constraint ----- %
% % Positions
% q_periodic = zeros(p.getNumPositions);
% q_periodic(2,2) = 1; %z
% q_periodic(3,3) = 1; %pitch-hip w/symmetry
% q_periodic(3,5) = 1; %pitch-hip w/symmetry
% q_periodic(4,6) = 1; %knee w/symmetry
% q_periodic(6,4) = 1; %knee w/symmetry
% q_periodic(5,5) = -1; %hip w/symmetry

% %Velocities
% v_periodic = zeros(p.getNumPositions);
% v_periodic(1,1) = 1; %x-vel
% v_periodic(2,2) = 1; %z-vel
% v_periodic(3,3) = 1; %pitch-hip w/symmetry
% v_periodic(3,5) = 1; %pitch-hip w/symmetry
% v_periodic(4,6) = 1; %knee w/symmetry
% v_periodic(6,4) = 1; %knee w/symmetry
% v_periodic(5,5) = -1; %hip w/symmetry
% 
% %Difference matrix to calculate velocities
% Diff = [-eye(p.getNumPositions), eye(p.getNumPositions)];

%Constraint Matrix
% R_periodic = [q_periodic, zeros(p.getNumPositions), zeros(p.getNumPositions), blkdiag(0, -eye(p.getNumPositions-1));
%               v_periodic*Diff, -eye(p.getNumPositions)*Diff];
          
% R_periodic = [q_periodic, blkdiag(0, -eye(p.getNumPositions-1))];
% 
% periodic_constraint = LinearConstraint(zeros(p.getNumPositions,1),zeros(p.getNumPositions,1),R_periodic);
% periodic_constraint = periodic_constraint.setName('periodicity');

% ----- Initial Guess ----- %

q0 = xtraj.eval(0);
qm = xtraj.eval(xtraj.tspan(2)/2);
qf = xtraj.eval(xtraj.tspan(2));

t_init = linspace(0,xtraj.tspan(2),N);
traj_init.x = xtraj; %PPTrajectory(foh([0 T0/2 T0],[q0 qm qf]));
traj_init.u = utraj; %PPTrajectory(zoh(t_init,randn(3,N)));
T_span = [1 T0];

traj_opt = VariationalTrajectoryOptimization(p,N,T_span);
traj_opt = traj_opt.addRunningCost(@running_cost_fun);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(q0(1:6)),1);
%traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qm(1:6)),8);
traj_opt = traj_opt.addStateConstraint(ConstantConstraint(qf(1:6)),N);
%traj_opt = traj_opt.addStateConstraint(periodic_constraint,{[1 N]});

% traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100);
% traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);
[xtraj,utraj,phitraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);

v = p.constructVisualizer();
v.playback(xtraj);

function [f,df] = running_cost_fun(h,x,u)
  f = h*u'*u;
  df = [u'*u zeros(1,12) 2*h*u'];
end

end
