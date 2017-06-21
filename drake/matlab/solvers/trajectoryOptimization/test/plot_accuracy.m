close all

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
sim_traj = PPTrajectory(foh(sim_traj.getBreaks(), sim_traj.eval(sim_traj.getBreaks())));

N = 2:2:20;
rms_pos1 = zeros(length(N),1);
rms_vel = zeros(length(N),1);

nq=plant.getNumPositions();

for i = 1:length(N)
  N(i)
  xtraj1 = variationalBrick(plant,N(i),x0);
  xtraj2 = contactImplicitBrick(plant,N(i),x0);
  %xtraj3 = contactImplicitBrick(plant,N(i),x0,'backward');
  ts1 = xtraj1.getBreaks();
  dt1 = ts1(2)-ts1(1);
  
%   xtraj = contactImplicitBrick(plant,N(i),x0);
%   xtraj_knots = xtraj.eval(xtraj.getBreaks());
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks()+(0.5*tf/N(i)),xtraj_knots)); % use foh vel
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks(),xtraj_knots)); % use foh vel

  x1 = xtraj1.eval(ts1);
  x2 = xtraj2.eval(ts1);
  %x3 = xtraj3.eval(ts1);
  q_true = sim_traj.eval(ts1);
  %v_true = sim_traj.eval(ts1+.5*dt1);
  rms_pos1(i) = rms(vec(x1(1:nq,:)-q_true(1:nq,:)));
  rms_pos2(i) = rms(vec(x2(1:nq,:)-q_true(1:nq,:)));
  %rms_pos3(i) = rms(vec(x3(1:nq,:)-q_true(1:nq,:)));
  %rms_vel(i) = rms(vec(x(nq+(1:nq),:)-v_true(nq+(1:nq),:)));
  figure(1)
  plot(N(1:i),rms_pos1(1:i),'bx-');
  hold on;
  plot(N(1:i),rms_pos2(1:i),'rx-');
  %plot(N(1:i),rms_pos3(1:i),'gx-');
  legend('Variational','Euler');
  
  figure(2);
  for j=1:6
    subplot(2,3,j);
    plot(ts1,x1(j,:),'b');
    hold on;
    plot(ts1,x2(j,:),'r');
    %plot(ts1,x3(j,:),'g');
    plot(ts1,q_true(j,:),'k--');
    hold off;
  end
  legend('Variational','Euler','True');
  
%   figure(3);
%   for j=1:6
%     subplot(2,3,j);
%     plot(ts,x(6+j,:),'b');
%     hold on;
%     plot(ts,v_true(6+j,:),'r');
%     hold off;
%   end
end

% figure();
% semilogy(N,rms_pos1);
% hold on
% semilogy(N,rms_pos2);

save('variational-rms.mat','N','rms_pos','rms_vel');
