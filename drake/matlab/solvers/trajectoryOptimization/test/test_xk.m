function [f,df,xdn,dfdx,dfdu,dfdSig] = test_xk(t,x,sigma,u,K,obj,pp)
% x = reshape(x,12,[]);
% sigma = Sig(1:obj.nx);
% u

%% obj
u_fdb_k = u - K*(sigma - x);

[xdn,df] = obj.plant.update(t,sigma,u_fdb_k);
f = sum(xdn);

%% gradient
[H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(sigma(1:obj.nx/2),sigma(obj.nx/2+1:obj.nx));
Hinv = inv(H);
Bmatrix = [1;zeros(5,1)];%B;hand coding

dfdu = [obj.plant.timestep^2*Hinv*Bmatrix;obj.plant.timestep*Hinv*Bmatrix];
dfdSig = df(:,2:obj.nx+1) - dfdu*K;
dfdx = dfdu*K;

switch pp
    case 1
        df = sum(dfdx,1)';
    case 2
        df = sum(dfdu(:));
    case 3
        df = sum(dfdSig,1)';
end

