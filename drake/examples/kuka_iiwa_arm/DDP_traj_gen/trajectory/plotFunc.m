function plotFunc

clear all

tmp = addpathTemporary(fullfile(pwd,'..'));

%% read arm data
%filename = 'joint_trajectory_interpolated.csv';

load joint_trajectory_interpolated.csv;
joint_pos = joint_trajectory_interpolated(:,1:7);
joint_vel = joint_trajectory_interpolated(:,8:14);

figure(1)
subplot(3,1,1)
plot(joint_vel(:,2))
subplot(3,1,2)
plot(joint_vel(:,4))
subplot(3,1,3)
plot(joint_vel(:,6))

a  =1
end