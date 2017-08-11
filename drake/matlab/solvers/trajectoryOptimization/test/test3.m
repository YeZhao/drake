function [f,dff] = test3(t,x,u_fdb_k,obj)
x = reshape(x,12,[]);

[xdn,df] = obj.plant.update(t,x,u_fdb_k);
[xdn1,df1] = obj.plant.update(t,xdn,u_fdb_k);
[xdn2,df2] = obj.plant.update(t,xdn1,u_fdb_k);
f = sum(xdn2);


dfdx = df(:,2:13);
% dfdx = sum(dfdx,1)';

dfdx1 = df1(:,2:13);
dfdx1 = sum(dfdx1,1);


dff = dfdx1*dfdx;
% dff = dff';
dff = dff(:);

