%=================================================
%  Date:        January 13th, 2016
%  File name:   friction_estimation.m
%  Author:      Ye Zhao (HarvardAgileRoboticsLab)
%=================================================

clear all
clc 
close all
fclose all
%% 
% sequence: p4, p7, p75, p8, p85, p9, 1, 
path = '~/kuka-dev/drake/drake/examples/kuka_iiwa_arm/experiment_data/friction_model/joint4/vel_p8';

read_joint_status_file;

% crop the beginning and ending part of data sequence
beginning_index = 50;
ending_index = length(joint_position_measured) - 1000;
joint_index = 4;

vertical_pos_index = [];
joint_vel_vertical_pos = [];
joint_torque_vertical_pos = [];

% collect all the data when the arm crosses zero vertical position
for i = beginning_index:ending_index
    if (sign(joint_position_measured(joint_index,i)) ~= sign(joint_position_measured(joint_index,i+1)))
        vertical_pos_index = [vertical_pos_index, i];
        joint_vel_vertical_pos = [joint_vel_vertical_pos, joint_velocity_measured(joint_index,i)];
        joint_torque_vertical_pos = [joint_torque_vertical_pos, joint_torque_measured(joint_index,i)];
    end
end

% delete the last element if the size of vertical_pos_index is odd
if mod(length(vertical_pos_index), 2) == 1
    vertical_pos_index = vertical_pos_index(1:end-1);
    joint_vel_vertical_pos = joint_vel_vertical_pos(1:end-1);
    joint_torque_vertical_pos = joint_torque_vertical_pos(1:end-1);
end

% classify position and negative sets for two directional motions.
joint_torque_positive_set = joint_torque_vertical_pos(joint_torque_vertical_pos>0);
joint_torque_negative_set = joint_torque_vertical_pos(joint_torque_vertical_pos<0);
joint_vel_positive_set = joint_vel_vertical_pos(joint_vel_vertical_pos>0);
joint_vel_negative_set = joint_vel_vertical_pos(joint_vel_vertical_pos<0);

if ((length(joint_torque_positive_set) ~= length(joint_torque_negative_set)) || (length(joint_vel_positive_set) ~= length(joint_vel_negative_set)))
    error('Error: sizes of posive and negative sets are not equal!');
end

% average the position and torque data
joint_torque_positive_avg = mean(joint_torque_positive_set);
joint_torque_negative_avg = mean(joint_torque_negative_set);
joint_vel_positive_avg = mean(joint_vel_positive_set);
joint_vel_negative_avg = mean(joint_vel_negative_set);

% For debugging
joint_torque_positive_avg
joint_torque_negative_avg
joint_vel_positive_avg
joint_vel_negative_avg
(joint_torque_positive_avg - joint_torque_negative_avg)/2
(joint_vel_positive_avg - joint_vel_negative_avg)/2

% friction_positive_set6 = [];
% friction_negative_set6 = [];
% 
% friction_positive_set6 = load('FRICTION_POSITIVE_SET6.dat'); 
% friction_negative_set6 = load('FRICTION_NEGATIVE_SET6.dat'); 
% 
% friction_positive_set6 = [friction_positive_set4;joint_vel_positive_avg,joint_torque_positive_avg];
% friction_negative_set6 = [friction_negative_set4;joint_vel_negative_avg,joint_torque_negative_avg];
% 
% save FRICTION_POSITIVE_SET6.dat friction_positive_set6 -ASCII
% save FRICTION_NEGATIVE_SET6.dat friction_negative_set6 -ASCII

% figure font 
bigTextSize = 20;
mediumTextSize = 18;
textSize = 16;
smallTextSize = 14;
legendMargin = 0.3;

% figure(1)
% coeffs_positive = polyfit(friction_positive_set6(:,1),friction_positive_set6(:,2), 1);
% x1 = linspace(0, 1.0, 100);
% y1 = polyval(coeffs_positive, x1);
% coeffs_negative = polyfit(friction_negative_set6(:,1),friction_negative_set6(:,2), 1);
% x2 = linspace(-1,0, 100);
% y2 = polyval(coeffs_negative, x2);
% plot(x1, y1,'Linewidth',2)
% hold on
% plot(x2, y2,'Linewidth',2)
% hold on
% plot(friction_positive_set6(:,1),friction_positive_set6(:,2),'ro','MarkerSize',10)
% hold on
% plot(friction_negative_set6(:,1),friction_negative_set6(:,2),'ro','MarkerSize',10)
% xlabel('joint velocity [rad/s]','fontsize',mediumTextSize)
% ylabel('joint torque [Nm]','fontsize',mediumTextSize)
% title('Joint 6 friction model','fontsize',bigTextSize)
% xlim([-1.2, 1.2])
% ylim([-2, 2])
% grid on
% %print -depsc Joint_6_friction_data
 
% figure(2)
% % data with manual shifting
% coeffs1 = polyfit(friction_positive_set6(:,1),friction_positive_set6(:,2), 1);
% x1 = linspace(0, 1.0, 100);
% y1 = polyval(coeffs1, x1);
% % x2 = linspace(-1,0, 100);
% % y2 = polyval(coeffs, x2);
% plot(x1, y1,'Linewidth',2)
% hold on
% % plot(x2, y2,'Linewidth',2)
% % hold on
% plot(friction_positive_set6(:,1),friction_positive_set6(:,2),'ro','MarkerSize',10)
% xlabel('joint velocity [rad/s]','fontsize',mediumTextSize)
% ylabel('joint torque [Nm]','fontsize',mediumTextSize)
% title('Joint 6 friction model','fontsize',bigTextSize)
% xlim([-1.2, 1.2])
% ylim([-2, 2])
% grid on
% %print -depsc Joint_6_friction_data

figure(3)
plot(joint_velocity_measured(joint_index,beginning_index:ending_index),joint_torque_measured(joint_index,beginning_index:ending_index))
grid on
xlabel('joint velocity [rad/s]')
ylabel('joint torque [Nm]')
title('Joint 6 Torque V.S. Velocity')
%print -depsc Joint_6_Torque_VS_Velocity

figure(4)
plot(joint_position_measured(joint_index,beginning_index:ending_index),joint_torque_measured(joint_index,beginning_index:ending_index))
grid on
xlabel('joint position [rad]')
ylabel('joint torque [Nm]')
title('Joint 6 Torque V.S. Position')
%print -depsc Joint_6_Torque_VS_POsition

figure(5)
plot(joint_position_measured(joint_index,beginning_index:ending_index),joint_velocity_measured(joint_index,beginning_index:ending_index))
grid on
xlabel('joint position [rad]')
ylabel('joint velocity [rad/s]')
title('Joint 6 Position V.S. Velocity')
%print -depsc Joint_6_Velocity_VS_Position

%joint 6
friction_positive_set6 = [0.1, 0.0840;
                          0.2, 0.0928;
                          0.3, 0.1025; 
                          0.5,    0.1213; %done
                          0.6,    0.1389;
                          0.65,   0.1350;
                          0.7,    0.1447;
                          0.75,   0.1514;
                          0.80,   0.1597;
                          0.85,   0.1651;
                          0.90,   0.1670;
                          1,      0.1800];

% saved data
% %joint 2
% friction_positive_set2 = [0.2039, 1.7140;
%                           0.3, 2.0362;
%                           0.4067, 2.5670;
%                           %0.4861, 3.7735;
%                           0.7478, 3.3125;
%                           0.8911, 8.3235;
%                           1.0411, 12.31;
%                           1.1885, 14.5450;
%                           1.3464, 14.8650;]

% %joint 4
% friction_positive_set4 = [0.1015, 0.4908;
%                           0.2018, 0.6111;
%                           0.3032, 0.8037; 
%                           0.4067, 0.8777; 
%                           %0.5030, 3.7735;
%                           0.6004, 1.2775;
%                           0.6494, 1.3065;
%                           %0.7079, 0.6660;
%                           0.7573, 1.5140;
%                           0.8088, 1.5415;
%                           0.8529, 1.6345;
%                           0.9001, 1.5795;]
%                           %1.0169, 0.3802;



% joint_vel_vertical_pos =
% 
%    -0.5018    0.5020   -0.5018    0.5020   -0.5017    0.5019   -0.5018    0.5019   -0.5019    0.5020   -0.5018    0.5019
% 
% joint_torque_vertical_pos
% 
% joint_torque_vertical_pos =
% 
%     0.0773   -0.2605    0.0707   -0.2523    0.0751   -0.2440    0.0689   -0.2442    0.0862   -0.2429    0.0844   -0.2610