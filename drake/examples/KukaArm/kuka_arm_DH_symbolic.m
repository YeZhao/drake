syms theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 theta_8 obj_x obj_y obj_z obj_yaw obj_pitch obj_roll T xB_x xB_y xB_z
state_vec = [theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 theta_8 obj_x obj_y obj_z obj_yaw obj_pitch obj_roll];

d1 = 0.4200;
d2 = 0.4000;
a = zeros(1,7);
alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
d = [0,0,d1,0,d2,0,0];
T = eye(4);

A(:,:,1) = [cos(theta_1), 0,sin(theta_1), 0;
    sin(theta_1),  0, -cos(theta_1), 0;
    0, 1, 0, 0;
    0, 0, 0, 1];

A(:,:,2) = [cos(theta_2), 0, -sin(theta_2), 0;
    sin(theta_2),  0, cos(theta_2), 0;
    0, -1, 0, 0;
    0, 0, 0, 1];

A(:,:,3) = [cos(theta_3), 0, -sin(theta_3), 0;
    sin(theta_3), 0, cos(theta_3), 0;
    0, -1, 0, 0.42;
    0, 0, 0, 1];

A(:,:,4) = [cos(theta_4), 0, sin(theta_4), 0;
    sin(theta_4), 0, -cos(theta_4), 0;
    0, 1, 0, 0;
    0, 0, 0, 1];
A(:,:,5) = [cos(theta_5), 0, sin(theta_5), 0;
    sin(theta_5), 0, -cos(theta_5), 0;
    0, 1, 0, 2/5;
    0, 0, 0,   1];

A(:,:,6) = [cos(theta_6), 0, -sin(theta_6), 0;
    sin(theta_6), 0, cos(theta_6), 0;
    0, -1, 0, 0;
    0, 0,  0, 1];

A(:,:,7) = [cos(theta_7), -sin(theta_7), 0, 0;
    sin(theta_7),  cos(theta_7), 0, 0;
    0,  0, 1, 0;
    0, 0, 0, 1];

for i=1:7
    T = T*A(:,:,i);
end
T = simplify(T);

T = [sin(theta_7)*(sin(theta_5)*(cos(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) - cos(theta_1)*sin(theta_2)*sin(theta_4)) - cos(theta_5)*(cos(theta_3)*sin(theta_1) + cos(theta_1)*cos(theta_2)*sin(theta_3))) - cos(theta_7)*(sin(theta_6)*(sin(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) + cos(theta_1)*cos(theta_4)*sin(theta_2)) + cos(theta_6)*(cos(theta_5)*(cos(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) - cos(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_3)*sin(theta_1) + cos(theta_1)*cos(theta_2)*sin(theta_3)))),   cos(theta_7)*(sin(theta_5)*(cos(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) - cos(theta_1)*sin(theta_2)*sin(theta_4)) - cos(theta_5)*(cos(theta_3)*sin(theta_1) + cos(theta_1)*cos(theta_2)*sin(theta_3))) + sin(theta_7)*(sin(theta_6)*(sin(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) + cos(theta_1)*cos(theta_4)*sin(theta_2)) + cos(theta_6)*(cos(theta_5)*(cos(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) - cos(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_3)*sin(theta_1) + cos(theta_1)*cos(theta_2)*sin(theta_3)))), sin(theta_6)*(cos(theta_5)*(cos(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) - cos(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_3)*sin(theta_1) + cos(theta_1)*cos(theta_2)*sin(theta_3))) - cos(theta_6)*(sin(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)) + cos(theta_1)*cos(theta_4)*sin(theta_2)), - (21*cos(theta_1)*sin(theta_2))/50 - (2*sin(theta_4)*(sin(theta_1)*sin(theta_3) - cos(theta_1)*cos(theta_2)*cos(theta_3)))/5 - (2*cos(theta_1)*cos(theta_4)*sin(theta_2))/5;
    cos(theta_7)*(sin(theta_6)*(sin(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) - cos(theta_4)*sin(theta_1)*sin(theta_2)) + cos(theta_6)*(cos(theta_5)*(cos(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) + sin(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_1)*cos(theta_3) - cos(theta_2)*sin(theta_1)*sin(theta_3)))) - sin(theta_7)*(sin(theta_5)*(cos(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) + sin(theta_1)*sin(theta_2)*sin(theta_4)) - cos(theta_5)*(cos(theta_1)*cos(theta_3) - cos(theta_2)*sin(theta_1)*sin(theta_3))), - cos(theta_7)*(sin(theta_5)*(cos(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) + sin(theta_1)*sin(theta_2)*sin(theta_4)) - cos(theta_5)*(cos(theta_1)*cos(theta_3) - cos(theta_2)*sin(theta_1)*sin(theta_3))) - sin(theta_7)*(sin(theta_6)*(sin(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) - cos(theta_4)*sin(theta_1)*sin(theta_2)) + cos(theta_6)*(cos(theta_5)*(cos(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) + sin(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_1)*cos(theta_3) - cos(theta_2)*sin(theta_1)*sin(theta_3)))), cos(theta_6)*(sin(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) - cos(theta_4)*sin(theta_1)*sin(theta_2)) - sin(theta_6)*(cos(theta_5)*(cos(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)) + sin(theta_1)*sin(theta_2)*sin(theta_4)) + sin(theta_5)*(cos(theta_1)*cos(theta_3) - cos(theta_2)*sin(theta_1)*sin(theta_3))),   (2*sin(theta_4)*(cos(theta_1)*sin(theta_3) + cos(theta_2)*cos(theta_3)*sin(theta_1)))/5 - (21*sin(theta_1)*sin(theta_2))/50 - (2*cos(theta_4)*sin(theta_1)*sin(theta_2))/5;
    sin(theta_7)*(sin(theta_5)*(cos(theta_2)*sin(theta_4) - cos(theta_3)*cos(theta_4)*sin(theta_2)) - cos(theta_5)*sin(theta_2)*sin(theta_3)) - cos(theta_7)*(cos(theta_6)*(cos(theta_5)*(cos(theta_2)*sin(theta_4) - cos(theta_3)*cos(theta_4)*sin(theta_2)) + sin(theta_2)*sin(theta_3)*sin(theta_5)) - sin(theta_6)*(cos(theta_2)*cos(theta_4) + cos(theta_3)*sin(theta_2)*sin(theta_4))),                                                                                                                                                                                                                                                                 cos(theta_7)*(sin(theta_5)*(cos(theta_2)*sin(theta_4) - cos(theta_3)*cos(theta_4)*sin(theta_2)) - cos(theta_5)*sin(theta_2)*sin(theta_3)) + sin(theta_7)*(cos(theta_6)*(cos(theta_5)*(cos(theta_2)*sin(theta_4) - cos(theta_3)*cos(theta_4)*sin(theta_2)) + sin(theta_2)*sin(theta_3)*sin(theta_5)) - sin(theta_6)*(cos(theta_2)*cos(theta_4) + cos(theta_3)*sin(theta_2)*sin(theta_4))),                                                                                                                                                            sin(theta_6)*(cos(theta_5)*(cos(theta_2)*sin(theta_4) - cos(theta_3)*cos(theta_4)*sin(theta_2)) + sin(theta_2)*sin(theta_3)*sin(theta_5)) + cos(theta_6)*(cos(theta_2)*cos(theta_4) + cos(theta_3)*sin(theta_2)*sin(theta_4)),                                                                        (21*cos(theta_2))/50 + (2*cos(theta_2)*cos(theta_4))/5 + (2*cos(theta_3)*sin(theta_2)*sin(theta_4))/5;
    0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                            1];


finger_contact_delta = 0.01;
right_finger_y_shift = 0.04;
R_DHbase_world = [-1,0,0;0,-1,0;0,0,1];% rotation from DH base coordinate to world coordinate (rotation around z axis)
pos_ee = R_DHbase_world*[T(1,4);T(2,4);T(3,4)+0.36];%0.36 is the base vertical position w.r.t. world coordinate
%R_ee = R_DHbase_world*T(1:3,1:3)*R_DHbase_world;
%R_ee_fwdkin = rpy2rotmat(iiwa_link_7_init(4:6));

fr1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
fr2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
fr3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
fr4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1025];

fl1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(theta_8-0.04);0.081+0.1225];
fl2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(theta_8-0.04);0.081+0.1025];
fl3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(theta_8-0.04);0.081+0.1225];
fl4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(theta_8-0.04);0.081+0.1025];

% Jacobian for A
JA_fr1_theta1 = diff(fr1,theta_1);
JA_fr1_theta2 = diff(fr1,theta_1);
JA_fr1_obj_yaw = diff(fr1,obj_yaw);

%dJ, second derivative
dJA_fr1_theta1_theta1 = diff(fr1,theta_1,2);
dJA_fr1_theta2_theta1 = diff(diff(fr1,theta_2),theta_1);
dJA_fr1_theta1_to_8_theta1 = diff([diff(fr1,theta_1);diff(fr1,theta_2);diff(fr1,theta_3);diff(fr1,theta_4);diff(fr1,theta_5);diff(fr1,theta_6);diff(fr1,theta_7);diff(fr1,theta_8)],theta_1);
dJA_fr2_theta1_to_8_theta1 = diff([diff(fr2,theta_1);diff(fr2,theta_2);diff(fr2,theta_3);diff(fr2,theta_4);diff(fr2,theta_5);diff(fr2,theta_6);diff(fr2,theta_7);diff(fr2,theta_8)],theta_1);

dJA_fr1_obj_yaw_theta1 = diff(diff(fr1,obj_yaw),theta_1);


b = [obj_x obj_y obj_z obj_yaw obj_pitch obj_roll]';

Rx = [1, 0, 0;0, cos(b(4)), -sin(b(4));0, sin(b(4)), cos(b(4))];
Ry = [cos(b(5)), 0, sin(b(5));0, 1, 0;-sin(b(5)), 0, cos(b(5))];
Rz = [cos(b(6)), -sin(b(6)), 0;sin(b(6)), cos(b(6)), 0;0, 0, 1];

R_world_to_B = Rz*Ry*Rx;

fl1 = R_world_to_B'*fl1;
fl2 = R_world_to_B'*fl2;
fl3 = R_world_to_B'*fl3;
fl4 = R_world_to_B'*fl4;

fr1 = R_world_to_B'*fr1;
fr2 = R_world_to_B'*fr2;
fr3 = R_world_to_B'*fr3;
fr4 = R_world_to_B'*fr4;

cylinder_radius = 0.03;
cylinder_height_half = 0.09;

b_local = R_world_to_B'*b(1:3);
contact_pt(1) = b_local(1) + cylinder_radius/norm(fl1(1:2)-b_local(1:2)) * (fl1(1) - b_local(1));
contact_pt(2) = b_local(2) + cylinder_radius/norm(fl1(1:2)-b_local(1:2)) * (fl1(2) - b_local(2));
contact_pt(3) = fl1(3);

% Jacobian for B
xA_ground = [cylinder_radius, -cylinder_radius, 0, 0;
    0, 0, cylinder_radius, -cylinder_radius;
    -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];

x_ground1 = b(1:3) + R_world_to_B*xA_ground(:,1);
x_ground2 = b(1:3) + R_world_to_B*xA_ground(:,2);
x_ground3 = b(1:3) + R_world_to_B*xA_ground(:,3);
x_ground4 = b(1:3) + R_world_to_B*xA_ground(:,4);

JB_ground1_theta1 = diff(x_ground1,theta_1);

right_normal1 = [fr1(1:2) - b_local(1:2);0];
right_normal1 = right_normal1./sqrt(right_normal1'*right_normal1);

% full version
fr1_B = cylinder_radius*right_normal1;
fr1_B(3) = fr1(3) - b_local(3);

%dJ, second derivative, full version
fr1_B = b(1:3) + R_world_to_B*fr1_B;
dJB_fr1_theta1 = diff(fr1_B,theta_1,2);
dJB_fr1_obj_yaw = diff(diff(fr1_B,obj_yaw),theta_1);

%simplified version, assuming the contact point is known
%fr1_B = b(1:3) + R_world_to_B*[xB_x xB_y xB_z]';
%JB_fr1_obj_yaw = diff(fr1_B,obj_yaw);

%dJ, second derivative, simplified version, assuming the contact point is known
fr1_B = b(1:3) + R_world_to_B*[xB_x xB_y xB_z]';
dJB_fr1_theta1_theta1 = diff(fr1_B,theta_1,2);
dJB_fr1_obj_yaw_theta1 = diff(diff(fr1_B,obj_yaw),theta_1);


