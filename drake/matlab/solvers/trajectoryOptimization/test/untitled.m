[xdn,df] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);

dfdt = df(:,1);
dfdx = df(:,2:13);
dfdu = df(:,14);

x0 = Sig(1:obj.nx,j,k);
x0 = x0(:);
[f1,dfdx1] = test(t,x0,u_fdb_k,obj);
ff = @(x) test(t,x,u_fdb_k,obj);
dd = DerivCheck(ff,x0)




%%
x0 = Sig(1:obj.nx,j,k);
x0 = x0(:);
[f1,dfdx1] = test2(t,x0,u_fdb_k,obj);
ff = @(x) test2(t,x,u_fdb_k,obj);

xxx0 = x0+randn(size(x0));
dd = DerivCheck(ff,xxx0)
% dd = DerivCheck(ff,x0)
norm(xxx0-x0)

%%
sigma0 = Sig(1:obj.nx,j,k);
u0 = 1;
x0 = Sig(1:obj.nx,j,k);

pp = 1;
[f,df,xdn] = test_xk(t,x0,sigma0,u0,K,obj,pp);
ff = @(x0) test_xk(t,x0,sigma0,u0,K,obj,pp);
dd = DerivCheck(ff,x0+randn(size(x0)))

pp = 2;
ff = @(u0) test_xk(t,x0,sigma0,u0,K,obj,pp);
dd = DerivCheck(ff,u0+randn(size(u0)))

pp = 3;
ff = @(sigma0) test_xk(t,x0,sigma0,u0,K,obj,pp);
dd = DerivCheck(ff,sigma0+randn(size(sigma0)))



%%
clc
% u1 = 1.7*10;
% x1 = randn(size(x1))*10;
% x0 = randn(size(x0))*10;
% u0 = randn*10;
% sigma0 = randn(size(sigma0))*10;

pp1 = 4;
[L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
ff = @(x1) test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
dd = DerivCheck(ff,x1+randn(size(x1)))


pp1 = 5;
[L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
ff = @(u1) test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
dd = DerivCheck(ff,u1+randn(size(u1)))

pp1 = 1;
[L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
ff = @(x0) test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
dd = DerivCheck(ff,x0+randn(size(x0)))


pp1 = 2;
[L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
ff = @(u0) test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
dd = DerivCheck(ff,u0+randn(size(u0)))

pp1 = 3;
[L,dL,sigma2] = test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
ff = @(sigma0) test_xk_1(t,x0,sigma0,u0,x1,u1,K,obj,pp1);
dd = DerivCheck(ff,sigma0+randn(size(sigma0)))






%%
clc
% u1 = 1.7*10;
% x1 = randn(size(x1))*10;
% x0 = randn(size(x0))*10;
% u0 = randn*10;
% sigma0 = randn(size(sigma0))*10;
% x2 = randn(size(x1));
% u2 = randn;

pp2 = 1;
[L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
ff = @(x0) test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
dd = DerivCheck(ff,x0+randn(size(x0)))

pp2 = 2;
[L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
ff = @(u0) test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
dd = DerivCheck(ff,u0+randn(size(u0)))

pp2 = 3;
[L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
ff = @(sigma0) test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
dd = DerivCheck(ff,sigma0+randn(size(sigma0)))

pp2 = 4;
[L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
ff = @(x1) test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
dd = DerivCheck(ff,x1+randn(size(x1)))


pp2 = 5;
[L,dL,sigma3] = test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
ff = @(u1) test_xk_2(t,x0,sigma0,u0,x1,u1,x2,u2,K,obj,pp2);
dd = DerivCheck(ff,u1+randn(size(u1)))




%%
clc
% N = 5;
% X0 = randn(N*12,1);
% U0 = randn(N,1);
% sigma_0 = randn(size(sigma0));

% ppn = 3;
% [L,dL,sigma_nn] = test_xk_n(t,X0,U0,sigma_0,K,obj,N,ppn);

ppn = 1;
ff = @(X0) test_xk_n(t,X0,U0,sigma_0,K,obj,N,ppn);
dd = DerivCheck(ff,X0+randn(size(X0)))

ppn = 2;
ff = @(U0) test_xk_n(t,X0,U0,sigma_0,K,obj,N,ppn);
dd = DerivCheck(ff,U0+randn(size(U0)))

ppn = 3;
ff = @(sigma_0) test_xk_n(t,X0,U0,sigma_0,K,obj,N,ppn);
dd = DerivCheck(ff,sigma_0)



% [ff,dff] = test_xk(t,X0(1:12),sigma_0,U0(1),K,obj,ppn);
% [ff,dff1,sigma3] = test_xk_1(t,X0(1:12),sigma_0,U0(1),X0(13:24),U0(2),K,obj,3);

% [ff,dff4,sigma3] = test_xk_1(t,X0(1:12),sigma_0,U0(1),X0(13:24),U0(2),K,obj,4);
% dff = [dff1;dff4];
% 
% [ff,dff,sigma3] = test_xk_2(t,X0(1:12),sigma_0,U0(1),X0(13:24),U0(2),X0(25:36),U0(3),K,obj,3);
% ff = @(sigma_0) test_xk_2(t,X0(1:12),sigma_0,U0(1),X0(13:24),U0(2),X0(25:36),U0(3),K,obj,3);
% dd = DerivCheck(ff,sigma0+randn(size(sigma0)))






















