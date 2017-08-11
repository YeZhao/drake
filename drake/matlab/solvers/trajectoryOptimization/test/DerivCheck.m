function dd = DerivCheck(funptr, X0, ~, varargin)

% DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
%`
%  Checks the analytic gradient of a function 'funptr' at a point X0, and
%  compares to numerical gradient.  Useful for checking gradients computed
%  for fminunc and fmincon.
%
%  Call with same arguments as you would call for optimization (fminunc).
%
% $id$

[~, JJ] = feval(funptr, X0, varargin{:});  % Evaluate function at X0
% rand('seed',1)

% Pick a random small vector in parameter space
tol = 1e-6;  % Size of numerical step to take
rr = randn(length(X0),1)*tol;  % Generate small random-direction vector

%                 rr_tmp = rr

%                 rr(13:24) = zeros(12,1);
%                 rr(1:36) = zeros(36,1);
%rr(49:end) = zeros(4,1);
%Derivs: Analytic vs. Finite Diff = [-4.204377749891e-02, -4.204377708857e-02]


% Evaluate at symmetric points around X0
f1 = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0
f2 = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0

% Print results
fprintf('Derivs: Analytic vs. Finite Diff = [%.4e, %.4e]\n', dot(rr, JJ), f2-f1);

dd =  dot(rr, JJ)-f2+f1;
end