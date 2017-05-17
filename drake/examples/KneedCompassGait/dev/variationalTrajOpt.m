function [p,xtraj,utraj,ltraj,ljltraj,z,F,info,traj_opt] = variationalTrajOpt()

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.ignore_self_collisions = true;
p = PlanarRigidBodyManipulator('../KneedCompassGait.urdf',options);

%todo: add joint limits, periodicity constraint

N = 10;
T = 3;
T0 = 3;

traj_opt = VariationalTrajectoryOptimization(p,N,T_span);

traj_opt = traj_opt.addRunningCost(@running_cost_fun);

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',100);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',200000);
[xtraj,utraj,ltraj,ljltraj,z,F,info] = traj_opt.solveTraj(t_init,traj_init);

v = p.constructVisualizer();
v.playback(xtraj);

function [f,df] = running_cost_fun(h,x,u)
  f = h*u'*u;
  df = [u'*u zeros(1,12) 2*h*u'];
end

end
