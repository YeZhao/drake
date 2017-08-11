function [L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1)
% x0 = reshape(x0,12,[]);
% sigma0 = Sig(1:obj.nx);
% u0

[f0,df0,sigma1,dfdx0,dfdu0,dfdSig0] = test_xk(t,x0,sigma0,u0,K,obj,pp1);
[f1,df1,sigma2,dfdx1,dfdu1,dfdSig1] = test_xk(t,x1,sigma1,u1,K,obj,pp1);
L = sum(sigma2);

%% obj
% u_fdb_k1 = u1 - K*(sigma1 - x1);
% [sigma2,df1] = obj.plant.update(t,sigma1,u_fdb_k1);
% L = sum(sigma2);
%
%
% %% gradient
% [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(sigma1(1:obj.nx/2),sigma1(obj.nx/2+1:obj.nx));
% Hinv = inv(H);
% Bmatrix = [1;zeros(5,1)];%B;hand coding
%
% dfdu1 = [obj.plant.timestep^2*Hinv*Bmatrix;obj.plant.timestep*Hinv*Bmatrix];
% dfdSig1 = df1(:,2:obj.nx+1) - dfdu1*K;
% dfdx1 = dfdu1*K;

%%
switch pp1
    case 1 % x0
        dL = sum(dfdSig1,1)*dfdx0;
        dL = dL(:);
    case 2 % u0
        dL = sum(dfdSig1,1)*dfdu0;
        dL = dL(:);
    case 3 % sigma0
        dL = sum(dfdSig1,1)*dfdSig0;
        dL = dL(:);
    case 4 % x1
        dL = sum(dfdx1,1)';
    case 5 % u1
        dL = sum(dfdu1(:));
end

