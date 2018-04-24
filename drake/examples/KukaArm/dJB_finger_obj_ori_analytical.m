function dJB_finger_obj_ori = dJB_finger_obj_ori_analytical(q,xB_x,xB_y,xB_z)
theta_1 = q(1);theta_2 = q(2);theta_3 = q(3);theta_4 = q(4);theta_5 = q(5);theta_6 = q(6);theta_7 = q(7);theta_8 = q(8);
obj_x = q(9); obj_y = q(10); obj_z = q(11); obj_yaw = q(12); obj_pitch = q(13); obj_roll = q(14);

% switch contact_index
%     case 1

dJB_finger_obj_ori = [ xB_y*(sin(obj_yaw)*sin(obj_roll) + cos(obj_yaw)*cos(obj_roll)*sin(obj_pitch)) + xB_z*(cos(obj_yaw)*sin(obj_roll) - cos(obj_roll)*sin(obj_yaw)*sin(obj_pitch)), cos(obj_yaw)*cos(obj_roll)*cos(obj_pitch)*xB_z - cos(obj_roll)*sin(obj_pitch)*xB_x + cos(obj_roll)*cos(obj_pitch)*sin(obj_yaw)*xB_y, xB_z*(cos(obj_roll)*sin(obj_yaw) - cos(obj_yaw)*sin(obj_roll)*sin(obj_pitch)) - xB_y*(cos(obj_yaw)*cos(obj_roll) + sin(obj_yaw)*sin(obj_roll)*sin(obj_pitch)) - cos(obj_pitch)*sin(obj_roll)*xB_x;
                        - xB_y*(cos(obj_roll)*sin(obj_yaw) - cos(obj_yaw)*sin(obj_roll)*sin(obj_pitch)) - xB_z*(cos(obj_yaw)*cos(obj_roll) + sin(obj_yaw)*sin(obj_roll)*sin(obj_pitch)), cos(obj_yaw)*cos(obj_pitch)*sin(obj_roll)*xB_z - sin(obj_roll)*sin(obj_pitch)*xB_x + cos(obj_pitch)*sin(obj_yaw)*sin(obj_roll)*xB_y, xB_z*(sin(obj_yaw)*sin(obj_roll) + cos(obj_yaw)*cos(obj_roll)*sin(obj_pitch)) - xB_y*(cos(obj_yaw)*sin(obj_roll) - cos(obj_roll)*sin(obj_yaw)*sin(obj_pitch)) + cos(obj_roll)*cos(obj_pitch)*xB_x;
                                                                                                                                 cos(obj_yaw)*cos(obj_pitch)*xB_y - cos(obj_pitch)*sin(obj_yaw)*xB_z,                                                           - cos(obj_pitch)*xB_x - cos(obj_yaw)*sin(obj_pitch)*xB_z - sin(obj_yaw)*sin(obj_pitch)*xB_y,                                                                                                                                                                                                                                                                                           0];

end