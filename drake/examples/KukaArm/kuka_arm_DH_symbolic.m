syms theta_1 theta_2 theta_3 theta_4 theta_5 theta_6 theta_7 T

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
 
