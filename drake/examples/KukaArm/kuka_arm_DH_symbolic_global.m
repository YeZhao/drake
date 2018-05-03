syms theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 theta_8 obj_x obj_y obj_z obj_yaw obj_pitch obj_roll T xB_x xB_y xB_z fr1x fr1y fr1z real
syms right_normal1_x right_normal1_y right_normal1_z real
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

J_fr1_theta_1 = jacobian(fr1,theta_1);
dJ_fr1_theta_1_theta_1 = jacobian(jacobian(fr1,theta_1),theta_1);


cylinder_radius = 0.03;
cylinder_height_half = 0.09;

fr1 = [fr1x;fr1y;fr1z];
b = [obj_x obj_y obj_z obj_yaw obj_pitch obj_roll]';
obj_pos = b(1:3);
obj_ori = b(4:6);
R_world_to_B = rpy2rotmat(obj_ori);

x_ground1A = obj_pos + R_world_to_B*[cylinder_radius;0;-cylinder_height_half];
x_ground2A = obj_pos + R_world_to_B*[-cylinder_radius;0;-cylinder_height_half];
x_ground3A = obj_pos + R_world_to_B*[0;cylinder_radius;-cylinder_height_half];
x_ground4A = obj_pos + R_world_to_B*[0;-cylinder_radius;-cylinder_height_half];

% Jacobian for A
JA_x_ground1_theta1 = diff(x_ground1A,theta_1);
JA_x_ground1_obj_x = diff(x_ground1A,obj_x);
JA_x_ground1_obj_yaw = diff(x_ground1A,obj_yaw);

JA_fr1_theta1 = diff(fr1,theta_1);
JA_fr1_theta2 = diff(fr1,theta_1);
JA_fr1_obj_yaw = diff(fr1,obj_yaw);

%dJ, second derivative
dJA_x_ground1_theta1 = diff([zeros(3*8,1);reshape(eye(3),[],1);diff(x_ground1A,obj_yaw);diff(x_ground1A,obj_pitch);diff(x_ground1A,obj_roll)],theta_1);

dJA_fr1_theta1_theta1 = diff(fr1,theta_1,2);
dJA_fr1_theta2_theta1 = diff(diff(fr1,theta_2),theta_1);
dJA_fr1_theta1_to_8_theta_1 = diff([diff(fr1,theta_1);diff(fr1,theta_2);diff(fr1,theta_3);diff(fr1,theta_4);diff(fr1,theta_5);diff(fr1,theta_6);diff(fr1,theta_7);diff(fr1,theta_8)],theta_1);
dJA_fr2_theta1_to_8_theta_1 = diff([diff(fr2,theta_1);diff(fr2,theta_2);diff(fr2,theta_3);diff(fr2,theta_4);diff(fr2,theta_5);diff(fr2,theta_6);diff(fr2,theta_7);diff(fr2,theta_8)],theta_1);
dJA_fr3_theta1_to_8_theta_1 = diff([diff(fr3,theta_1);diff(fr3,theta_2);diff(fr3,theta_3);diff(fr3,theta_4);diff(fr3,theta_5);diff(fr3,theta_6);diff(fr3,theta_7);diff(fr3,theta_8)],theta_1);
dJA_fr4_theta1_to_8_theta_1 = diff([diff(fr4,theta_1);diff(fr4,theta_2);diff(fr4,theta_3);diff(fr4,theta_4);diff(fr4,theta_5);diff(fr4,theta_6);diff(fr4,theta_7);diff(fr4,theta_8)],theta_1);
dJA_fl1_theta1_to_8_theta_1 = diff([diff(fl1,theta_1);diff(fl1,theta_2);diff(fl1,theta_3);diff(fl1,theta_4);diff(fl1,theta_5);diff(fl1,theta_6);diff(fl1,theta_7);diff(fl1,theta_8)],theta_1);
dJA_fl2_theta1_to_8_theta_1 = diff([diff(fl2,theta_1);diff(fl2,theta_2);diff(fl2,theta_3);diff(fl2,theta_4);diff(fl2,theta_5);diff(fl2,theta_6);diff(fl2,theta_7);diff(fl2,theta_8)],theta_1);
dJA_fl3_theta1_to_8_theta_1 = diff([diff(fl3,theta_1);diff(fl3,theta_2);diff(fl3,theta_3);diff(fl3,theta_4);diff(fl3,theta_5);diff(fl3,theta_6);diff(fl3,theta_7);diff(fl3,theta_8)],theta_1);
dJA_fl4_theta1_to_8_theta_1 = diff([diff(fl4,theta_1);diff(fl4,theta_2);diff(fl4,theta_3);diff(fl4,theta_4);diff(fl4,theta_5);diff(fl4,theta_6);diff(fl4,theta_7);diff(fl4,theta_8)],theta_1);

dJA_fr1_theta1_to_8_obj_x = diff([diff(fr1,obj_x);diff(fr1,obj_y);diff(fr1,obj_z);diff(fr1,obj_yaw);diff(fr1,obj_pitch);diff(fr1,obj_roll)],theta_3);


dJA_fr1_obj_yaw_theta1 = diff(diff(fr1,obj_yaw),theta_1);
dJA_fr1_obj_x_theta1 = diff(diff(fr1,obj_x),theta_1);
dJA_fr1_obj_y_theta1 = diff(diff(fr1,obj_y),theta_1);
dJA_fr1_obj_z_theta1 = diff(diff(fr1,obj_z),theta_1);

fl1 = R_world_to_B'*fl1;
fl2 = R_world_to_B'*fl2;
fl3 = R_world_to_B'*fl3;
fl4 = R_world_to_B'*fl4;

fr1 = R_world_to_B'*fr1;
fr2 = R_world_to_B'*fr2;
fr3 = R_world_to_B'*fr3;
fr4 = R_world_to_B'*fr4;

b_local = R_world_to_B'*b(1:3);
contact_pt(1,:) = b_local(1) + cylinder_radius/norm(fr1(1:2)-b_local(1:2)) * (fr1(1) - b_local(1));
contact_pt(2,:) = b_local(2) + cylinder_radius/norm(fr1(1:2)-b_local(1:2)) * (fr1(2) - b_local(2));
contact_pt(3,:) = fr1(3);
contact_pt = R_world_to_B*contact_pt;

% Jacobian for B
xB_ground = [cylinder_radius, -cylinder_radius, 0, 0;
    0, 0, cylinder_radius, -cylinder_radius;
    -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];

x_ground1B = b(1:3) + R_world_to_B*xB_ground(:,1);
x_ground2B = b(1:3) + R_world_to_B*xB_ground(:,2);
x_ground3B = b(1:3) + R_world_to_B*xB_ground(:,3);
x_ground4B = b(1:3) + R_world_to_B*xB_ground(:,4);

JB_ground1_theta1 = diff(x_ground1B,theta_1);

right_normal1 = [fr1(1:2) - b_local(1:2);0];
right_normal1 = right_normal1./sqrt(right_normal1'*right_normal1);

dright_normal1_dfr1 = jacobian(right_normal1,[fr1x,fr1y,fr1z]);
ddright_normal1_ddfr1x = jacobian(right_normal1,[fr1x]);
ddright_normal1_ddfr1y = jacobian(right_normal1,[fr1y]);
ddright_normal1_ddfr1z = jacobian(right_normal1,[fr1z]);

% chain rule
right_normal1 = [right_normal1_x right_normal1_y right_normal1_z]';

fr1_B = cylinder_radius*right_normal1;
fr1_B(3) = fr1(3) - b_local(3);

fr1_B = b(1:3) + R_world_to_B*fr1_B;
JB_fr1_theta1 = diff(fr1_B,theta_1);

% chain rule
dfr1_B_dright_normal = jacobian(fr1_B,[right_normal1_x right_normal1_y right_normal1_z]);
% dfr1_B_dright_normal = [ (3*cos(obj_roll)*cos(obj_pitch))/100, (3*cos(obj_roll)*sin(obj_yaw)*sin(obj_pitch))/100 - (3*cos(obj_yaw)*sin(obj_roll))/100, 0;
%  (3*cos(obj_pitch)*sin(obj_roll))/100, (3*cos(obj_yaw)*cos(obj_roll))/100 + (3*sin(obj_yaw)*sin(obj_roll)*sin(obj_pitch))/100, 0;
%               -(3*sin(obj_pitch))/100,                                                    (3*cos(obj_pitch)*sin(obj_yaw))/100, 0];
 

dvecAdx = jacobian(vec(dfr1_B_dright_normal),[fr1x;fr1y;fr1z]);

ddfr1_B_ddright_normal1 = jacobian(dfr1_B_dright_normal(:,1),[right_normal1_x right_normal1_y right_normal1_z]);
ddfr1_B_ddright_normal2 = jacobian(dfr1_B_dright_normal(:,2),[right_normal1_x right_normal1_y right_normal1_z]);
ddfr1_B_ddright_normal3 = jacobian(dfr1_B_dright_normal(:,3),[right_normal1_x right_normal1_y right_normal1_z]);

%[ddfr1_B_ddright_normal1*dright_normal1_dfr1(:,1),
%ddfr1_B_ddright_normal2*dright_normal1_dfr1(:,1),
%ddfr1_B_ddright_normal3*dright_normal1_dfr1(:,1)] * dright_normal1_dfr1

%+ dright_normal1_dfr1*ddright_normal1_ddfr1x
JB_fr1_B_fr1x = jacobian(fr1_B,[fr1x]);

JB_fr1_B_fr1orignial = jacobian(fr1_B,[fr1x,fr1y,fr1z]);
dJB_fr1_B_fr1orignial1 = jacobian(JB_fr1_B_fr1orignial(:,1),[fr1x,fr1y,fr1z]);
dJB_fr1_B_fr1orignial2 = jacobian(JB_fr1_B_fr1orignial(:,2),[fr1x,fr1y,fr1z]);
dJB_fr1_B_fr1orignial3 = jacobian(JB_fr1_B_fr1orignial(:,3),[fr1x,fr1y,fr1z]);

dJB_fr1_theta1 = diff(fr1_B,theta_1,2);
dJB_fr1_obj_yaw = diff(diff(fr1_B,obj_yaw),theta_1);

%simplified version, assuming the contact point is known
fr1_B = b(1:3) + R_world_to_B*[xB_x xB_y xB_z]';
JB_fr1_obj_yaw = diff(fr1_B,obj_yaw);
JB_fr1_obj_ori = [diff(fr1_B,obj_yaw),diff(fr1_B,obj_pitch),diff(fr1_B,obj_roll)];

%dJ, second derivative, simplified version, assuming the contact point is known
fr1_B = b(1:3) + R_world_to_B*[xB_x xB_y xB_z]';
dJB_fr1_theta1_theta1 = diff(fr1_B,theta_1,2);
dJB_fr1_obj_yaw_theta1 = diff(diff(fr1_B,obj_yaw),theta_1);
