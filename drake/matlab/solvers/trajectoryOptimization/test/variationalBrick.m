function variationalBrick()

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
plant = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
x0 = [0;0;.8;0.05*randn(3,1);zeros(6,1)];

N=6;
tf=.5;

prog = VariationalTrajectoryOptimization(plant,N,tf);

prog = prog.addStateConstraint(ConstantConstraint(x0(1:6)),1);
prog = prog.addStateConstraint(ConstantConstraint(x0(1:6)),2);

traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));

[xtraj,utraj,ltraj,z,F,info,infeasible] = prog.solveTraj(tf,traj_init);

v = constructVisualizer(plant);
v.playback(xtraj);

end

