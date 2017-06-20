function [xtraj,utraj,ctraj,btraj,straj,z] = variationalBrick(plant,N)

if nargin<1
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
%plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%x0 = [0; 1; .1; 0; 0; 0];
file = fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf');
plant = RigidBodyManipulator(file,options);
end
if nargin < 2
  N=10;
end

q0 = [0
      0
      1.0000
      0.1
     -0.1
      0.0];
    
x0 = [q0;0*q0];

tf=1.0;

t_init = linspace(0,tf,N);
traj_init.x = PPTrajectory(foh([0 tf],[x0, x0]));

options.s_weight = 100;

traj_opt = VariationalTrajectoryOptimization(plant,N,tf,options);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(q0),1);  
traj_opt = traj_opt.addVelocityConstraint(ConstantConstraint(0*q0),1);

traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',1000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = traj_opt.solveTraj(t_init,traj_init);
toc

v = constructVisualizer(plant);
v.playback(xtraj);

% ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

% sim_traj = ts_plant.simulate([0,tf],x0);
% 
% 
% ts = sim_traj.getBreaks();
% %tt = xtraj.getBreaks();
% xx = xtraj.eval(ts);
% xs = sim_traj.eval(ts);
% 
% figure(1);
% for i=1:12
%   subplot(2,6,i);
%   plot(ts,xx(i,:),'b');
%   hold on;
%   plot(ts,xs(i,:),'r');
%   hold off;
% end

end

