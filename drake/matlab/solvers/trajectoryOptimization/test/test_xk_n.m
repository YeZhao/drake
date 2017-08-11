function [L,dL,sigma_nn] = test_xk_n(t,X0,U0,sigma_0,K,obj,N,ppn)
% x0 = reshape(x0,12,[]);
% sigma0 = Sig(1:obj.nx);
% u0

% N = obj.N;
dL_u = ones(1,12,N);
dL_x = ones(12,12,N);
dL_sigma = 1;
sigma_nn = sigma_0;

for nn = 1:N
    
    %% hat(X_N+1)
    x_nn = X0((nn-1)*12+1:nn*12,:);
    u_nn = U0(nn);
    [~,~,sigma_nn_1,dfdx_nn,dfdu_nn,dfdSig_nn] = test_xk(t,x_nn,sigma_nn,u_nn,K,obj,ppn);
    
    dL_u(:,:,nn) = dfdu_nn';
    dL_x(:,:,nn) = dfdx_nn';
    dL_sigma = dL_sigma*dfdSig_nn';
    
    if nn>1
        uu = dL_u(:,:,1:nn-1);
        uu = permute(uu,[3,1,2]);
        [nx,ny,nz] = size(uu);
        uu = reshape(uu,[],size(uu,3));
        uu = uu*dfdSig_nn';
        uu = permute(reshape(uu,[nx,ny,nz]),[2,3,1]);
        dL_u(:,:,1:nn-1) = uu;
        
        xx = dL_x(:,:,1:nn-1);
        xx = permute(xx,[3,1,2]);
        [nx,ny,nz] = size(xx);
        xx = reshape(xx,[],size(xx,3));
        xx = xx*dfdSig_nn';
        xx = reshape(xx,[nx,ny,nz]);
        xx = permute(reshape(xx,[nx,ny,nz]),[2,3,1]);
        dL_x(:,:,1:nn-1) = xx;
        
    end
    
    %% X_N+1
    sigma_nn = sigma_nn_1;
end

dL_u = sum(permute(dL_u,[3,1,2]),3);
dL_x = sum(permute(dL_x,[3,1,2]),3);

L = sum(sigma_nn);

%%
switch ppn
    case 1 % X0
        dL = dL_x';
        dL = dL(:);
    case 2 % U0
        dL = dL_u;
    case 3 % sigma0
        dL = sum(dL_sigma,2);
end

