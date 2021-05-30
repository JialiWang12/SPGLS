function obj = fval(X,y,z,gamma,w)
m = size(X, 1);
X = [X, ones(m, 1)];
alpha = w' * w;
obj = norm(alpha/gamma * z + X * w - y - alpha/gamma * y)^2 / (1+alpha/gamma)^2;
