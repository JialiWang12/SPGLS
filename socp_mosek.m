function [w_opt, optval,time, timeeig,lambda] = socp_mosek(X,y,z,gamma,varargin)

[m, n] = size(X);
Y = [X ones(m,1)];
A11 = Y'*Y; A11=(A11'+A11)/2;
A12 = Y'*(z-y);
A13 = -Y'*y;
A22 = (norm(z-y))^2;
A23 = -(z-y)'*y;
A33 = y'*y;

A12_mid = 1/sqrt(gamma)*(A12-A13);
A11_bar = [A11, A12_mid ; A12_mid', 1/gamma *(A22+A33-2*A23)];
A12_bar = [A12+A13; 1/sqrt(gamma) * (A22-A33)];
A22_bar = A22+A33+2*A23;

A = [A11 A12 A13; A12' A22 A23; A13' A23' A33];
B = zeros(n+3); B(n+2:n+3, n+2:n+3) = ones(2,2);
C = zeros(n+3); C(1:n+1,1:n+1) = 1/gamma * speye(n+1); C(n+2,n+3) = -0.5; C(n+3,n+2) = -0.5; 

teig = tic;
[H,D] = eig(full(A11_bar));
timeeig = toc(teig);
A12_tilde = H' * A12_bar;
d = diag(D); 
vec_ones = ones(n+2,1);
idx_temp = reshape(3:3*(n+2)+2,n+2,3);
cone_idx = reshape(idx_temp',1,[]);
start_idx = 1:3:length(cone_idx);
temp0 = tic;
clear prob;
[r, res] = mosekopt('symbcon');
prob.c = [1 0 zeros(1,3*(n+2))];
c1 = [4, 1, zeros(1,n+2), ones(1,n+2),zeros(1,n+2)];
c2 = zeros(n+2,2+3*(n+2));c2(:,2) = -1/gamma;c2(:,3:n+4) =diag(ones(n+2,1));
c = [c1;c2];
prob.a = [c];
prob.buc = [A22_bar,d'];
prob.blc = [-inf,d'];
prob.blx = [-inf, -inf, zeros(1,n+2),-inf * ones(1,n+2),A12_tilde'.*sqrt(2) ];
prob.bux = [inf * ones(1,2+2*(n+2)),A12_tilde'.*sqrt(2)];
prob.cones.type = [repmat(res.symbcon.MSK_CT_RQUAD,n+2,1)];
prob.cones.sub = cone_idx;
prob.cones.subptr = start_idx;
temp1 = tic;
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;
[r,res] = mosekopt('maximize info echo(0)',prob,param);
% [r,res] = mosekopt('maximize info echo(0)');
ttotal = toc(temp0);time = ttotal;
tsolver = toc(temp1);

sol = res.sol.itr.xx';
mu = sol(1); lambda = sol(2);
optval = mu;


D_quad = A - mu*B + lambda*C; 
wful = [D_quad(1:n+3,1:n+3);[zeros(1,n+2) 1]]\[zeros(n+3,1); 1];
w_opt = wful(1:n+1);


