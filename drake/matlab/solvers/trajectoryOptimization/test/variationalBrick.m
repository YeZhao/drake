function [xtraj,utraj,ctraj,btraj,straj,z] = variationalBrick()

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
%plant = PlanarRigidBodyManipulator(fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf'),options);
%x0 = [0; 1; .1; 0; 0; 0];
file = fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf');
plant = RigidBodyManipulator(file,options);
x0 = [0;0;1;.1;0;0;zeros(6,1)];

N=11;
tf=1;

prog = VariationalTrajectoryOptimization(plant,N,tf);

%prog = prog.addStateConstraint(ConstantConstraint(x0(1:3)),1);
%prog = prog.addStateConstraint(ConstantConstraint(x0(1:3)+[0 -0.5*9.81*.1*.1 0]'),2);

prog = prog.addPositionConstraint(ConstantConstraint(x0(1:6)),1);
%prog = prog.addVelocityConstraint(ConstantConstraint(x0(7:end)),1);
prog = prog.addPositionConstraint(ConstantConstraint(x0(1:6)+[0 0 -0.5*9.81*.1*.1 0 0 0]'),2);

prog = prog.setSolverOptions('snopt','IterationsLimit',1000000);
% prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance', 1e-3);
% prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-3);
% prog = prog.setConstraintErrTol(1e-4);

traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));

tic
[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = prog.solveTraj(tf,traj_init);
toc

v = constructVisualizer(plant);
v.playback(xtraj);

end

