function w_opt = ridge_regression(X, y, nu)
    m = size(X, 1);
    n = size(X, 2);
    
    cvx_begin quiet
        variable w(n)
        variable b
        minimise ( sum_square_abs(X*w + b - y) + nu*w.'*w ) 
    cvx_end
    
    w_opt = [w; b];
end