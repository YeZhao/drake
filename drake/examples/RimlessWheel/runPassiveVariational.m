function runPassiveVariational

options.floating = true;
options.terrain = RigidBodyFlatTerrain();
w = warning('off','Drake:RigidBody:SimplifiedCollisionGeometry');
p = PlanarRigidBodyManipulator('RimlessWheel.urdf',options);
warning(w);
q0 = [0;1+rand;randn];
v0 = [5*rand;randn;5*rand];
x0 = p.resolveConstraints([q0;v0]);
q0 = x0(1:3);
q1 = q0 + .1*v0;
x1 = p.resolveConstraints([q1; v0]);
q1 = x1(1:3);

N = 21;
Tf = 2;

trajOpt = VariationalTrajectoryOptimization(p,N,[Tf Tf]);
trajOpt = trajOpt.addStateConstraint(ConstantConstraint(q0),1);
trajOpt = trajOpt.addStateConstraint(ConstantConstraint(q1),2);
trajOpt = trajOpt.setSolverOptions('snopt','IterationsLimit',200000);

t_init = 0:.1:2;
traj_init.x = PPTrajectory(foh([0 Tf],[q0 q1]));

[xtraj,utraj,ctraj,btraj,straj,z,F,info,infeasible_constraint_name] = trajOpt.solveTraj(t_init,traj_init);

v = p.constructVisualizer();
v.axis = [0 5 -.1 3];
v.playback(xtraj);
