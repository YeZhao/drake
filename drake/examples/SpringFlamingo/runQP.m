
options.use_bullet = false;
options.twoD = true;
options.view = 'right';
options.floating = true;
options.terrain = RigidBodyFlatTerrain();
s = 'urdf/spring_flamingo_passive_ankle.urdf';
dt = 0.001;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
r = TimeSteppingRigidBodyManipulator(s,dt,options);
warning(w);

v = r.constructVisualizer;

data=load('data/spring_2_steps.mat');
data.xtraj = data.xtraj.setOutputFrame(r.getStateFrame);
v.playback(data.xtraj);
x0=data.xtraj.eval(0);
c = IDQP(r,data.xtraj,data.utraj,data.ctraj);

sys = feedback(r,c);
% Forward simulate dynamics with visulazation, then playback at realtime
S=warning('off','Drake:DrakeSystem:UnsupportedSampleTime');
output_select(1).system=1;
output_select(1).output=1;
sys = mimoCascade(sys,v,[],[],output_select);
warning(S);
xtraj = simulate(sys,[0 data.xtraj.tspan(2)],x0);
playback(v,xtraj,struct('slider',true));


keyboard