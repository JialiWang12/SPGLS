function [w_opt, optval, time,lambda] = singlesdp_mosek(X,y, z, gamma, varargin)
[m,n] = size(X);
Y = [X ones(m,1)];
A11 = Y'*Y; A11=(A11'+A11)/2;
A12 = Y'*(z-y);
A13 = -Y'*y;
A22 = (norm(z-y))^2;
A23 = -(z-y)'*y;
A33 = y'*y;
A = [A11 A12 A13; A12' A22 A23; A13' A23' A33];
B = zeros(n+3); B(n+2:n+3, n+2:n+3) = ones(2,2);
C = zeros(n+3); C(1:n+1,1:n+1) = 1/gamma*speye(n+1); C(n+2,n+3) = -0.5; C(n+3,n+2) = -0.5; 


% begin mosek
temp0 = tic;
clear prob
[r, res] = mosekopt('symbcon');

% coefficients in the objective
[r0 c0 v0] = find(tril(A));
% coefficients in the constrains
[r1 c1 v1] = find(tril(B));
[r2 c2 v2] = find(tril(C));
% idx in the objective
prob.bardim = [n+3];
prob.barc.subj = ones(length(v0),1);
prob.barc.subk = [r0];
prob.barc.subl = [c0];
prob.barc.val = [v0];
% idx in the constrain
prob.blc = [1,0];
prob.buc = [1,0];
prob.a = sparse([], [], [], 2, 0);
prob.bara.subi = [ones(length(v1),1);2 * ones(length(v2),1)];
prob.bara.subj = [ones(length(v1)+length(v2),1)];
prob.bara.subk = [r1;r2];
prob.bara.subl = [c1;c2];
prob.bara.val = [v1;v2];

param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;
[r,res] = mosekopt('minimize info echo(0)',prob, param);
ttotal = toc(temp0); time = ttotal;
optval = res.info.MSK_DINF_INTPNT_PRIMAL_OBJ;
mu = res.sol.itr.y(1);
lambda = -res.sol.itr.y(2);

D = A - mu*B + lambda*C;
wful = [D(1:n+3,1:n+3);[zeros(1,n+2) 1]]\[zeros(n+3,1); 1];
w_opt = wful(1:n+1);
