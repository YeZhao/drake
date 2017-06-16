% a = [0,0,1]';
% b = [0,1,0]';

n = 3;%dimension
v = randn(1,n);
v = v./sqrt(v*v');
a = v';
a = [0,0,1]';

v = randn(1,n);
v = v./sqrt(v*v');
b = v';

v = cross(a,b);

s = norm(v);
c = dot(a,b);

    v_cross_prod = [0, -v(3), v(2);
                v(3), 0, -v(1);
                -v(2), v(1), 0];
            
R = eye(3) + v_cross_prod + v_cross_prod*v_cross_prod*(1-c)/s^2;

%R = (1 -)

b = R*a;


syms n1 n2 n3

D1 = [(n2^2+n1^2*n3)/(n2^2+n1^2)];

dD1 = gradient(D1, [n1,n2,n3]);
