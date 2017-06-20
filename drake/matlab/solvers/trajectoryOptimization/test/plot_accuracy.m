

options.terrain = RigidBodyFlatTerrain();
options.floating = true;

file = fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf');
plant = RigidBodyManipulator(file,options);

q0 = [0
      0
      1.5000
      0
      0
      0];
    
v0 = zeros(6,1);
v0(1) = 10.0;
x0 = [q0;v0];

tf=1.0;

ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

sim_traj = ts_plant.simulate([0,tf],x0);
ts = sim_traj.getBreaks();
xs = sim_traj.eval(ts);

N = 2:2:30;
rms_pos = zeros(length(N),1);
rms_vel = zeros(length(N),1);

nq=plant.getNumPositions();

for i = 1:length(N)
  N(i)
  xtraj = variationalBrick(plant,N(i),x0);
%   xtraj = contactImplicitBrick(plant,N(i),x0);
  xtraj_knots = xtraj.eval(xtraj.getBreaks());
  xtraj2 = PPTrajectory(foh(xtraj.getBreaks()+(0.5*tf/N(i)),xtraj_knots)); % use foh vel
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks(),xtraj_knots)); % use foh vel
  xp = xtraj.eval(ts);
  xv = xtraj2.eval(ts);
  rms_pos(i) = rms(vec(xp(1:nq,:)-xs(1:nq,:)));
  rms_vel(i) = rms(vec(xv(nq+(1:nq),:)-xs(nq+(1:nq),:)));
  figure(1)
  subplot(1,2,1);
  plot(N(1:i),rms_pos(1:i),'b');
  subplot(1,2,2);
  plot(N(1:i),rms_vel(1:i),'r');
  
  figure(2);
  for j=1:6
    subplot(2,3,j);
    plot(ts,xv(6+j,:),'b');
    hold on;
    plot(ts,xs(6+j,:),'r');
    hold off;
  end
end

save('variational-rms.mat','N','rms_pos','rms_vel');
