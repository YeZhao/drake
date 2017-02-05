function paramEstSyntheticData_KukaArm_J5J6J7_dynamic_test
clear all

tmp = addpathTemporary(fullfile(pwd,'..'));

%% read arm data
%path = '~/kuka-dev/drake/drake/examples/kuka_iiwa_arm/experiment_data/friction_model/joint5_6_7_dynamic_motion_nice_data';
path = '~/kuka-dev-estimation/drake/drake/examples/kuka_iiwa_arm/experiment_data/dynamic_test/J4_chirp';

read_joint_status_file;
joint_index = 6;

% beginning_index_set = [];
% ending_index_set = [];
% beginning_index_trigger = 1;
% ending_index_trigger = 0;
% 
% % extract the trajectory segment in static states
% for i = 1:length(joint_velocity_measured_raw)
%     if abs(joint_velocity_measured_raw(joint_index,i)) < 1e-6 && beginning_index_trigger == 1
%         beginning_index_set = [beginning_index_set,i];
%         beginning_index_trigger = 0;
%         ending_index_trigger = 1;
%     end
%     if abs(joint_velocity_measured_raw(joint_index,i)) > 1e-6 && ending_index_trigger == 1
%         if i - beginning_index_set(end) < 50
%             beginning_index_set = beginning_index_set(1:end-1);% ignore the bad set, sometimes the velocity becomes smaller than the threthold but then crosses it again. 
%         else
%             ending_index_set = [ending_index_set,i];
%         end
%         ending_index_trigger = 0;
%         beginning_index_trigger = 1;
%     end
% end
% 
% if length(beginning_index_set) > length(ending_index_set)
%     beginning_index_set = beginning_index_set(1:end-1);% crop the last element of beginning index set
% end
% 
% % % crop the beginning and ending part of data sequence
% % beginning_index = 1000;%2800;%1000;%200;
% % ending_index = 1800;%length(joint_position_measured) - 400;%1800;%1800;%3600;%2000;%
% 
joint_begin_index = 5;
joint_end_index = 7;
% 
% joint_velocity_measured = [];
% joint_position_measured = [];
% joint_acceleration_measured = [];
% joint_torque_measured = [];
% 
% beginning_index_set = beginning_index_set(1:10);
% ending_index_set = ending_index_set(1:10);
% 
% for i=1:length(beginning_index_set)
%     joint_velocity_measured = [joint_velocity_measured, joint_velocity_measured_raw(joint_begin_index:joint_end_index,beginning_index_set(i):ending_index_set(i))];
%     joint_position_measured = [joint_position_measured, joint_position_measured_raw(joint_begin_index:joint_end_index,beginning_index_set(i):ending_index_set(i))];
%     joint_acceleration_measured = [joint_acceleration_measured, joint_acceleration_measured_raw(joint_begin_index:joint_end_index,beginning_index_set(i):ending_index_set(i))];
%     joint_torque_measured = [joint_torque_measured, joint_torque_measured_raw(joint_begin_index:joint_end_index,beginning_index_set(i):ending_index_set(i))];
% end

beginning_index = 1;
ending_index = length(joint_position_measured_raw);

joint_velocity_measured = joint_velocity_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_position_measured = joint_position_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_acceleration_measured = joint_acceleration_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = joint_torque_measured_raw(joint_begin_index:joint_end_index,beginning_index:ending_index);
joint_torque_measured = -joint_torque_measured;

span = 0.005;
for i=1:joint_end_index-joint_begin_index+1
    joint_acceleration_measured(i,:) = smooth(joint_acceleration_measured(i,:),span,'loess');
    joint_velocity_measured(i,:) = smooth(joint_velocity_measured(i,:),span,'loess');
    joint_torque_measured(i,:) = smooth(joint_torque_measured(i,:),span,'loess');
end

%%
% Synthetic parameter estimation mode
% 'base':       Base case parameter estimation - intial params = true params
% 'paramerr':   w/ parameter error but no measurement noise
% 'measnoise':  w/ parameter error and w/ measurement noise
% 'delay':      w/ parameter error and w/ measurement noise and delay (Not complete)
mode = 'measnoise';

% Parameter Estimation model
% 'dynamic'     = use dynamic model - requires qdd
% 'energetic'   = use energetic model - doesn't require qdd
parameterEstimationOptions.model = 'dynamic';

% Parameter Estimation robot
% 'Acrobot'
% 'KukaArm'
parameterEstimationOptions.robot = 'KukaArm';

% Method by which to obtain qdd (not used if using energetic model)
% 'manipul':   Use acrobot manipulator equations to estimate true qdd
% 'derivative': Take the derivative of qd
qddmode = 'manipul';

% Parameter Estimation method
% 'nonlinprog'  = nonlinear least squares (LS) to solve problem
% 'linprog'       = linear LS on lumped params then nonlinear LS to recover
%                 original parameters
% 'simerr'      = minimize simulation error
% 'lsqnonlin'   = use MATLAB's built-in nonlinear least squares solver (debugging)
parameterEstimationOptions.method = 'nonlinprog';

% Option to print from estimator
% 'noprint'     = Do not print output from parameterEstimation.m, but will 
%                 still print output of paramEstSyntheticData.m script
% 'printEst'    = only print estimated from parameterEstimation.m
% 'printAll'	= print estimated and original from parameterEstimation.m
parameterEstimationOptions.print_result = 'noprint';

% Standard deviation of the data input error
if strcmp(parameterEstimationOptions.model,'energetic')
    % In this case, for theta1,theta2,theta1dot,theta2dot
    noisestd = sqrt([.0005, .0005, .0005, .0005, .0005, .0005, .0005, .0007, .0007, .0007, .0007, .0007, .0007, .0007]);
else
    % In this case, for theta1,theta2,theta1dot,theta2dot,theta1doubledot,theta2doubledot
    noisestd = sqrt([.0005, .0005, .0005, .0005, .0005, .0005, .0005, .0007, .0007, .0007, .0007, .0007, .0007, .0007, .0012, .0012, .0012, .0012, .0012, .0012, .0012]);
end
% Standard deviation of the parameter value percent error
paramstd = 1/5;

r = KukaArmPlant_J5J6J7_dynamic_test;
rtrue = KukaArmPlant_J5J6J7_dynamic_test;
    
% if ~strcmp(mode,'base')
%     % Perturb original parameter estimates with random percentage error
%     % normally distributed with standard dev = paramstd, and greater than -1
%     paramerr = randn(1,10)*paramstd;
%     while sum(paramerr<=-1)~=0
%         paramerr(paramerr<-1) = randn(1,sum(paramerr<-1))*paramstd;
%     end
% 
%     % rtrue.l1 = rtrue.l1 + rtrue.l1*paramerr(1); 
%     % rtrue.l2 = rtrue.l2 + rtrue.l2*paramerr(2); 
%     % rtrue.m1 = rtrue.m1 + rtrue.m1*paramerr(3); 
%     % rtrue.m2 = rtrue.m2 + rtrue.m2*paramerr(4);
%     rtrue.bv1_positive  = rtrue.bv1_positive + rtrue.bv1_positive*paramerr(5);
%     rtrue.bv2_positive  = rtrue.bv2_positive + rtrue.bv2_positive*paramerr(6);
%     rtrue.c1x = rtrue.c1x + rtrue.c1x*paramerr(7); 
%     rtrue.c2x = rtrue.c2x + rtrue.c2x*paramerr(8); 
%     rtrue.I1xx = rtrue.I1xx + rtrue.I1xx*paramerr(9);  
%     rtrue.I2xx = rtrue.I2xx + rtrue.I2xx*paramerr(10);
% end

outputFrameNames = r.getOutputFrame.getCoordinateNames();

% %% Test on swingup up data
% % [utraj,xtraj] = swingUpTrajectory(rtrue);
% % Ts = .01; breaks=getBreaks(utraj); T0 = breaks(1); Tf = breaks(end);
% % tsamples = T0:Ts:Tf;
% Ts = .01; 
% tsamples = 0:Ts:6;
% utraj = zeros(601,7);
% xtraj = zeros(601,14);
% xsamples = zeros(601,14);%eval(xtraj,tsamples)';
% usamples = zeros(601,7);%eval(utraj,tsamples)';

% % calcuate qddot
% nq = r.num_positions;
% qdd = zeros(length(tsamples),nq);
% for i=1:length(tsamples)
%     [H,C,B] = manipulatorDynamics(rtrue,xsamples(i,1:nq)',xsamples(i,nq+(1:nq))');
%     qdd(i,:) = (H\(B*usamples(i,:)' - C))';
% end

% %% Generate second derivative
% if ~strcmp(parameterEstimationOptions.model,'energetic')
%     nq = r.num_positions;
%     if strcmp(qddmode,'manipul')
%         qdd = zeros(length(tsamples),nq);
%         for i=1:length(tsamples)
%           [H,C,B] = manipulatorDynamics(rtrue,xsamples(i,1:nq)',xsamples(i,nq+(1:nq))');
%           qdd(i,:) = (H\(B*usamples(i,:)' - C))';
%         end
%     else
%         % Differentiating to get the second derivative of the state variables
%         % TODO: Try lowpass filter
%         qdd = deriv(xtraj,T0:Ts:Tf)'; qdd = qdd(:,nq+(1:nq));
%     end
%     xsamples = [xsamples qdd];
%     outputFrameNames = [outputFrameNames;'theta1doubledot';'theta2doubledot'];
%     
% %     qddTrue = qdd;
% %     qdd = deriv(xtraj,T0:Ts:Tf)'; qdd = qdd(:,nq+(1:nq));
% %     figure;plot(qddTrue(:,1),'-g'); hold on; plot(qdd(:,1),'-r')
% %     title(['q_1 True Acceleration vs. Derivative of Velocity']);
% %     xlabel('Sample Time')
% %     ylabel({'Acceleration Magnitude' '(m/s)/sample'})
% %     legend('True Acceleration','Derivative of Vel');
% end

% Only for simulated data
% %% Add gaussian noise to measurements
% if strcmp(mode,'measnoise') || strcmp(mode,'delay')
%     measurementNoise = randn(size(xsamples))*diag(noisestd);
% else
%     measurementNoise = 0;
% end
% xsamplesfinal = xsamples+measurementNoise;

arm_state = [joint_position_measured',joint_velocity_measured',joint_acceleration_measured'];
arm_input = joint_torque_measured';

% % Debugging purposes
% v = AcrobotVisualizer(r);
% vtrue = AcrobotVisualizer(rtrue);
% vtrue.playback(xtraj);
% 
% % plot(xsamples(:,1)); hold on; plot(xsamplesfinal(:,1));
% plot(xsamples); hold on; plot(xsamplesfinal);
% % Testing to see if estimate is produced by local minima
% r = rtrue;

outputFrameNames = [outputFrameNames;'theta5doubledot';'theta6doubledot';'theta7doubledot'];
Ts = .01;
data = iddata(arm_state,arm_input,Ts,'InputName',r.getInputFrame.getCoordinateNames(),'OutputName',outputFrameNames);
[estimated_parameters, estimated_delay] = parameterEstimation(r,data,parameterEstimationOptions);

%% Print out results
coords = getCoordinateNames(r.getParamFrame);
p_true = double(rtrue.getParams);
p_init = double(r.getParams);
fprintf('\nParameter estimation results:\n\n');
fprintf('  Param  \tKuka CAD Model    \tEstimated\n');
fprintf('  -----  \t--------\t\t---------\n');
for i=1:length(coords)
  fprintf('%7s  \t%8.6f\t\t%8.6f\n',coords{i},p_true(i),estimated_parameters(i));
end
