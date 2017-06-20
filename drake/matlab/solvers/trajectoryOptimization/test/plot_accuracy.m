

options.terrain = RigidBodyFlatTerrain();
options.floating = true;

file = fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf');
plant = RigidBodyManipulator(file,options);

q0 = [0
      0
      1.0000
      0.1
     -0.1
      0.0];
    
x0 = [q0;0*q0];

tf=1.0;

ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

sim_traj = ts_plant.simulate([0,tf],x0);
ts = sim_traj.getBreaks();
xs = sim_traj.eval(ts);

N = 3:20;
y = zeros(length(N),1);
for i = 1:length(N)
  N(i)
%   xtraj = variationalBrick(plant,N(i));
  xtraj = contactImplicitBrick(plant,N(i));
  xtraj_knots = xtraj.eval(xtraj.getBreaks());
  xtraj2 = PPTrajectory(foh(xtraj.getBreaks(),xtraj_knots)); % use foh vel
  xx = xtraj2.eval(ts);
  y(i) = rms(vec(xx-xs));
  figure(1)
  plot(N(1:i),y(1:i));
  
end

save('data-rms.mat','N','y');
