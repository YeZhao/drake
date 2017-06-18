q = [0,0.25,0.8573,0.45]';

%sqrt(1 - 0.3^2 - 0.2^2)
w = q(1);x = q(2);y = q(3);z = q(4);

R = [1-2*y^2-2*z^2, 2*x*y-2*z*w, 2*x*z+2*y*w;
     2*x*y+2*z*w, 1-2*x^2-2*z^2, 2*y*z-2*x*w;
     2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x^2-2*y^2];
 
n = R(:,3);

mu_q = [0,0.5,0.7071,0.5]';
mu_w = mu_q(1);mu_x = mu_q(2);mu_y = mu_q(3);mu_z = mu_q(4);

Rbar = [1-2*mu_y^2-2*mu_z^2, 2*mu_x*mu_y-2*mu_z*mu_w, 2*mu_x*mu_z+2*mu_y*mu_w;
        2*mu_x*mu_y+2*mu_z*mu_w, 1-2*mu_x^2-2*mu_z^2, 2*mu_y*mu_z-2*mu_x*mu_w;
        2*mu_x*mu_z-2*mu_y*mu_w, 2*mu_y*mu_z+2*mu_x*mu_w, 1-2*mu_x^2-2*mu_y^2];
 
F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
     -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
     2*mu_w, -2*mu_x, -2*mu_y, 2*mu_z];

% F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
%      -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
%      0, -4*mu_x, -4*mu_y, 0];

Fc = Rbar(:,3) - F*mu_q;

nbar_r = F*q + Fc;

% Quaternion to Euler angle conversion
psi = atan2(2*(q(1)*q(2) + q(3)*q(4)),1-2*(q(2)^2 + q(3)^2));
theta = asin(2*(q(1)*q(3) - q(4)*q(2)));
phi = atan2(2*(q(1)*q(4) + q(2)*q(3)), 1-2*(q(3)^2 + q(4)^2));

q_Euler = [psi,theta,phi]';

R_Euler = [cos(phi)*cos(theta), -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi), sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
           sin(phi)*cos(theta), cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), -cos(phi)*sin(psi)+ sin(phi)*sin(theta)*cos(psi);
           -sin(theta), cos(theta)*sin(psi), cos(theta)*cos(psi)];

% Quaternion to Euler angle conversion
mu_psi = atan2(2*(mu_q(1)*mu_q(2) + mu_q(3)*mu_q(4)),1-2*(mu_q(2)^2 + mu_q(3)^2));
mu_theta = asin(2*(mu_q(1)*mu_q(3) - mu_q(4)*mu_q(2)));
mu_phi = atan2(2*(mu_q(1)*mu_q(4) + mu_q(2)*mu_q(3)), 1-2*(mu_q(3)^2 + mu_q(4)^2));

mu_q_Euler = [mu_psi,mu_theta,mu_phi]';

Rbar_Euler = [cos(mu_phi)*cos(mu_theta), -sin(mu_phi)*cos(mu_psi) + cos(mu_phi)*sin(mu_theta)*sin(mu_psi), sin(mu_phi)*sin(mu_psi)+cos(mu_phi)*sin(mu_theta)*cos(mu_psi);
           sin(mu_phi)*cos(mu_theta), cos(mu_phi)*cos(mu_psi) + sin(mu_phi)*sin(mu_theta)*sin(mu_psi), -cos(mu_phi)*sin(mu_psi)+ sin(mu_phi)*sin(mu_theta)*cos(mu_psi);
           -sin(mu_theta), cos(mu_theta)*sin(mu_psi), cos(mu_theta)*cos(mu_psi)];

%Euler angle formulation
F_Euler = [sin(phi)*cos(psi) - cos(phi)*sin(theta)*sin(psi), cos(phi)*cos(theta)*cos(psi), cos(phi)*sin(psi) - sin(phi)*sin(theta)*cos(psi);
           -cos(phi)*cos(psi) - sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(theta)*cos(psi), sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
           -cos(theta)*sin(psi), -sin(theta)*cos(psi), 0];

Fc_Euler = Rbar_Euler(:,3) - F_Euler*mu_q_Euler;

nbar_r_Euler = F_Euler*q_Euler + Fc_Euler;


disp('the end')