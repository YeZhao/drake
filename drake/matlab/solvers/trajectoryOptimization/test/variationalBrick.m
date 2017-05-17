function [xtraj,utraj,phitraj,ctraj,btraj,straj,z] = variationalBrick()

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
plant = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%plant = RigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickNoContact.urdf'),options);
x0 = [0;0;1;0*randn(3,1);zeros(6,1)];

N=7;
tf=.6;

prog = VariationalTrajectoryOptimization(plant,N,tf);

prog = prog.addStateConstraint(ConstantConstraint(x0(1:6)),1);
prog = prog.addStateConstraint(ConstantConstraint(x0(1:6)+[0 0 -0.5*9.81*.1*.1 0 0 0]'),2);

prog = prog.setSolverOptions('snopt','IterationsLimit',1000000);
% prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance', 1e-3);
% prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-3);
% prog = prog.setConstraintErrTol(1e-4);

traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));

[xtraj,utraj,phitraj,ctraj,btraj,straj,z,F,info,infeasible] = prog.solveTraj(tf,traj_init);
v = constructVisualizer(plant);
v.playback(xtraj);

end

