q = [0,0.3,0.2,0.9327]';

%sqrt(1 - 0.3^2 - 0.2^2)
w = q(1);x = q(2);y = q(3);z = q(4);

R = [1-2*y^2-2*z^2, 2*x*y-2*z*w, 2*x*z+2*y*w;
     2*x*y+2*z*w, 1-2*x^2-2*z^2, 2*y*z-2*x*w;
     2*x*z-2*y*w, 2*y*z+2*x*w, 1-2*x^2-2*y^2];
 
n = R(:,3);

mu_q = [0,0,0,1]';
mu_w = mu_q(1);mu_x = mu_q(2);mu_y = mu_q(3);mu_z = mu_q(4);

Rbar = [1-2*mu_y^2-2*mu_z^2, 2*mu_x*mu_y-2*mu_z*mu_w, 2*mu_x*mu_z+2*mu_y*mu_w;
     2*mu_x*mu_y+2*mu_z*mu_w, 1-2*mu_x^2-2*mu_z^2, 2*mu_y*mu_z-2*mu_x*mu_w;
     2*mu_x*mu_z-2*mu_y*mu_w, 2*mu_y*mu_z+2*mu_x*mu_w, 1-2*mu_x^2-2*mu_y^2];
 
F = [2*mu_y, 2*mu_z, 2*mu_w, 2*mu_x;
     -2*mu_x, -2*mu_w, 2*mu_z, 2*mu_y;
     2*mu_w, -2*mu_x, -2*mu_y, 2*mu_z];
 
Fc = Rbar(:,3) - F*mu_q;

nbar_r = F*q + Fc;

disp('the end')