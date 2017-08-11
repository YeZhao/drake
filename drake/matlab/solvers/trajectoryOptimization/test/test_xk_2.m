function [L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2)
% x0 = reshape(x0,12,[]);
% sigma0 = Sig(1:obj.nx);
% u0

[f0,df0,sigma1,dfdx0,dfdu0,dfdSig0] = test_xk(t,x0,sigma0,u0,K,obj,pp2);
[f1,df1,sigma2,dfdx1,dfdu1,dfdSig1] = test_xk(t,x1,sigma1,u1,K,obj,pp2);
[f2,df2,sigma3,dfdx2,dfdu2,dfdSig2] = test_xk(t,x2,sigma2,u2,K,obj,pp2);
L = sum(sigma3);

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
switch pp2
    case 1 % x0
        dL = sum(dfdSig2*dfdSig1*dfdx0,1);
        dL = dL(:);
    case 2 % u0
        dL = sum(dfdSig2*dfdSig1*dfdu0,1);
        dL = dL(:);
    case 3 % sigma0
        dL = sum(dfdSig2*dfdSig1*dfdSig0,1);
        dL = dL(:);
    case 4 % x1
        dL = sum(dfdSig2*dfdx1,1);
    case 5 % u1
        dL = sum(dfdSig2*dfdu1,1);
    case 6 % x2
        dL = sum(dfdx2,1)';
    case 7 % u2
        dL = sum(dfdu2(:));
end

