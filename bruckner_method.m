function [w,iter] = bruckner_method(X, y, z, gamma)
    m = size(X, 1);
    n = size(X, 2);
    obj = bruckner_objective([X, ones(m, 1)], y, z, gamma);
    nonlincon = bruckner_constraints([X, ones(m, 1)], y, z, gamma);
    v0 = zeros(n+1+m, 1);
    hessian_finder = get_hessian_finder([X, ones(m, 1)], y, z, gamma);
    options = optimoptions('fmincon','Algorithm','interior-point', 'MaxFunctionEvaluations', 1e+10, 'Display', 'off', 'SpecifyConstraintGradient', false,'SpecifyObjectiveGradient', false, 'OptimalityTolerance', 1e-2);
    [x_opt,fval,exitflag,output] = fmincon(obj, v0, [], [], [], [], [], [], nonlincon, options);
    w = x_opt(1:n+1);
    iter = output.iterations;
end

function obj = bruckner_objective(X, y, z, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    function [loss, g] = bruckner(v)
        w = reshape(v(1:n), n, 1);
        tau = v(n+1: n+m);
        tau = reshape(tau, m, 1);
        inner = X*w + w.'*w.*tau - y;
        loss = norm(inner, 2)^2;
        
        if nargout  > 1
            inner = X*w + (w.'*w)*tau - y;
            gw = 2*X.'*inner + 2*inner.'*tau*w + 2*tau.'*inner*w;
            gt = 2*(w.'*w)*inner;
            g = [gw; gt];
        end
    end
    
    obj = @bruckner;
end

function nonlincon = bruckner_constraints(X, y, z, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    function [c, ceq, gc, gceq] = bruckner(v)
        w = reshape(v(1:n), n , 1);
        tau = v(n+1: n+m);
        tau = reshape(tau, m, 1);
        inner_term = X*w + w.'*w.*tau;
        ceq = 2*(m/gamma)*(inner_term - z) + tau;
        c = [];
        if nargout > 2
            gc = [];
            gceq_w = 2*(m /gamma)*X + 4*(m/gamma)*tau*w.';
            gceq_tau = eye(m) + 2*(m/gamma)*w.'*w*eye(m);
            gceq = [gceq_w, gceq_tau].';
        end
    end

    nonlincon = @bruckner;
end



function map_values = heatmap(gamma, q_start)
    X_tr = randn(20, 2);
    y_tr = X_tr*[0.5; 0.5] + 1.0;
    X_te = randn(50, 2);
    y_te = X_te*[0.5; 0.5] + 1.0;
    
    w1_list = linspace(-5, 5, 32);
    w2_list = linspace(-5, 5, 32);
    
    simSpace = [32, 32];
    numSims = prod(simSpace);
    bruckner_values = zeros(numSims, 1);
    bisect_values = zeros(numSims, 1);
    
    
    for idx=1:numSims
        [i, j] = ind2sub(simSpace, idx);
        z_tr = X_tr*[w1_list(i); w2_list(j)] + 1.0;
        z_te = X_te*[w1_list(i); w2_list(j)] + 1.0;
        w_b = bruckner_method(X_tr, y_tr, z_tr, gamma);
        bruckner_loss = compute_loss(X_te, y_te, z_te, gamma, w_b);
        bruckner_values(idx)= bruckner_loss;
        fprintf('In heatmap bruck %d \n', idx);
    end
    
    for idx=1:numSims
        [i, j] = ind2sub(simSpace, idx);
        z_tr = X_tr*[w1_list(i); w2_list(j)] + 1.0;
        z_te = X_te*[w1_list(i); w2_list(j)] + 1.0;
        w_s = bisect(X_tr,y_tr, z_tr, gamma, 10.0, q_start, 0.001, fid ,solver);
        bisect_loss = compute_loss(X_te, y_te, z_te, gamma, w_s);
        bisect_values(idx)= bisect_loss;
        fprintf('In heatmap bisect %d \n', idx);
    end
    
    map_values = bruckner_values - bisect_values;
end

function map_values = heatmap_ridge(gamma, q_start, fid, solver)
    X_tr = randn(20, 2);
    y_tr = X_tr*[0.5; 0.5] + 1.0;
    X_te = randn(50, 2);
    y_te = X_te*[0.5; 0.5] + 1.0;
    
    w1_list = linspace(-5, 5, 32);
    w2_list = linspace(-5, 5, 32);
    
    simSpace = [32, 32];
    numSims = prod(simSpace);
    bisect_values = zeros(numSims, 1);
    
    
    for idx=1:numSims
        [i, j] = ind2sub(simSpace, idx);
        z_tr = X_tr*[w1_list(i); w2_list(j)] + 1.0;
        z_te = X_te*[w1_list(i); w2_list(j)] + 1.0;
        w_s = bisect(X_tr,y_tr, z_tr, gamma, 10.0, q_start, 0.001, fid, solver);
        bisect_values(idx) = compute_loss(X_te, y_te, z_te, gamma, w_s);
        fprintf('In heatmap bisect %d \n', idx);
    end
    
    simSpace = [32, 32];
    numSims = prod(simSpace);
    ridge_values = zeros(numSims, 8);
    nu_range = [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100 ,1000];
    
    for idx=1:numSims
        [i, j] = ind2sub(simSpace, idx);
        z_tr = X_tr*[w1_list(i); w2_list(j)] + 1.0;
        z_te = X_te*[w1_list(i); w2_list(j)] + 1.0;
        for j=1:8
            nu = nu_range(j);
            w_r = ridge_regression(X_tr, y_tr, nu);
            ridge_loss = compute_loss(X_te, y_te, z_te, gamma, w_r);
            ridge_values(idx, j)= ridge_loss;
            fprintf('In heatmap ridge %d \n', idx);
        end
    end
    
    ridge_values = min(ridge_values, 2);
    
    map_values = ridge_values - bisect_values;
end

function hessian_obj = get_obj_hessian(X, y, z, gamma)
    m = size(X, 1);
    n = size(X, 2);
    function hessian = hessian_obj_inner(v)
        w = v(1:n);
        tau = v(n+1:m+n);
        dwdw = 2*X.'*X + 4*X.'*tau*w.' + 4*w*tau.'*X + 8*tau.'*tau*w*w.' + 4*(X*w + w.'*w*tau - y).'*tau*eye(n);
        dtdt = 2*(w.'*w)*(w.'*w)*eye(m);
        dwdt = 2*X'*(w.'*w) + 4*(w.'*w)*w*tau.' + 4*w*(X*w + (w.'*w)*tau - y).';
        
        hessian = [ dwdw, dwdt; dwdt.', dtdt];
    end

    hessian_obj = @hessian_obj_inner;
end

function hessian_ci = get_hessian_cons(X, y, z, gamma)
    m = size(X, 1);
    n = size(X, 2);
    function hessian_c = hessian_ci_inner(v, i)
        w = v(1:n);
        tau = v(n+1:m+n);
        t = tau(i);
        dwdw = 4*(m / gamma)*t*eye(n);
        dwdt = zeros(n, m);
        dwdt(:, i) = 2*(m / gamma)*w;
        dtdt = zeros(m, m);
        
        hessian_c = [dwdw, dwdt; dwdt.', dtdt];
    end
    hessian_ci = @hessian_ci_inner;
end

function hessian_finder = get_hessian_finder(X, y, z, gamma)
    hessian_obj = get_obj_hessian(X, y, z, gamma);
    hessian_ci = get_hessian_cons(X, y, z, gamma);
    m = size(X, 1);
    n = size(X, 2);
    
    function hessian = find_hessian(v, lambda)
        hessian = hessian_obj(v);
        le = lambda.eqnonlin;
        for i=1:m
            hessian = hessian + le(i)*hessian_ci(v, i);
        end
    end

    hessian_finder = @find_hessian;
    
end

