x_full = x_full +randn(size(x_full))*0.1;
Fext_full = Fext_full +randn(size(Fext_full))*0.1;

X0 = [x_full; Fext_full];
x0 = x_full;
fun = @(x) robustVariancecost_check_xonly(obj, x, X0)
dd0 = DerivCheck(fun, x0)



ii = 1;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd1 = DerivCheck(fun, x_ii)

ii = 2;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd2 = DerivCheck(fun, x_ii)

ii = 3;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd3 = DerivCheck(fun, x_ii)


dd = [dd0;dd1;dd2;dd3];



ii = 4;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd4 = DerivCheck(fun, x_ii)

dd = [dd0;dd1;dd2;dd3;dd4];


ii = 5;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd5 = DerivCheck(fun, x_ii)


dd = [dd0;dd1;dd2;dd3;dd4;dd5];









%%

ii = 4;
fun = @(x) robustVariancecost_check_xonly_ii(obj, x, X0, ii);
x_ii = X0((ii-1)*12+1:ii*12);
x_ii = x_ii+randn(size(x_ii))*0.1;
dd4 = DerivCheck(fun, x_ii)
[c,dc]=fun(x_ii);
% c_numeric = c;
% dc_numeric = dc;
[c_numeric,dc_numeric] = geval(fun,x_ii,struct('grad_method','numerical'));

valuecheck(dc',dc_numeric,1e-5);
valuecheck(c,c_numeric,1e-5);








