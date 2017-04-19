function fallingBrickLCP

options.floating = 'quat';
options.terrain = RigidBodyFlatTerrain();
% options.ignore_self_collisions = true;
% options.use_bullet = false;
s = 'FallingBrickContactPoints.urdf';
% s = 'FallingBrickBetterCollisionGeometry.urdf';
p = TimeSteppingRigidBodyManipulator_full_with_constraint_rho5_new(s,.01,options);
p = p.addRobotFromURDF(s,[],[],options);
% x0 = [0;1;2;rpy2quat(randn(3,1));randn(6,1)];
fix1 = [0.4400; 0.1017; 2.7873];
fix2 = [-1.1667; -1.8543; -1.1407];
fix3 = [-0.2185    0.5413    0.3893    0.7512    1.7783    1.2231   -1.2833   -2.3290    0.9019   -1.8356    0.0668    0.0355]';
% x0 = [0;1;2;rpy2quat(randn(3,1));2;1;2;rpy2quat(randn(3,1));randn(12,1)];
x0 = [0;1;2;rpy2quat(fix1);2;1;2;rpy2quat(fix2);fix3];
x0 = p.resolveConstraints(x0);

if 0 
  v = p.constructVisualizer();
  sys = cascade(p,v);
  sys.simulate([0 8],x0);
  return;
end

v = p.constructVisualizer();
v.drawWrapper(0,x0);
xtraj = p.simulate([0 4],x0);
v.playback(xtraj);