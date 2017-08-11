function [f,dfdx] = test(t,x,u_fdb_k,obj)
x = reshape(x,12,[]);
[xdn,df] = obj.plant.update(t,x,u_fdb_k);
f = sum(xdn);


dfdx = df(:,2:13);
dfdx = sum(dfdx,1)';
