function [bisect_means, bruckner_means, ridge_best_means, ssdp_means, socp_means] = get_mse(X, y, z, q_start, bias_upper, tol, gamma_range, nu_range, const, varargin)
    if nargin == 10
         bruckner_flag = varargin{1};
    end
    m = size(X, 1);
    num_gamma = length(gamma_range);
    num_nu = length(nu_range);
    rng(2021);
    CVO = cvpartition(m, 'Kfold', 10);
    
    % initialization
    bisect_means = zeros(num_gamma,1); bruckner_means = zeros(num_gamma,1);
    ridge_best_means = zeros(num_gamma,1); ssdp_means = zeros(num_gamma,1); socp_means = zeros(num_gamma,1);
    % flatten out
    simSpace = [10, num_gamma, num_nu];
    numSims = prod(simSpace);
    ridge_results = zeros(numSims, 1);

    % normalization
    X = normalize(X,'range');
    nor = const*max(abs(y));
    y = y/nor;
    z = z/nor;
    
    for i=1:10
        trIdx = CVO.training(i);
        teIdx = CVO.test(i);
        X_tr = X(trIdx, :);
        y_tr = y(trIdx);
        z_tr = z(trIdx);
        X_te = X(teIdx, :);
        y_te = y(teIdx);
        z_te = z(teIdx);

        n = size(X_tr, 2);

        % ridge regression method        
        fprintf('Doing ridge regression\n');
        ridge_params = ridge(y_tr, X_tr, nu_range, 0);
        % Need to swap first and last rows of ridge params as constant is
        % put at top of matrix
        row_1 = ridge_params(1, :);
        row_n1 = ridge_params(n+1, :);
        ridge_params(1, :) = row_n1;
        ridge_params(n+1,:) = row_1;
        for j=1:num_gamma
            gamma = gamma_range(j);
            for k=1:num_nu
                wr = ridge_params(:, k);
                % reverse normalization loss
                loss  = compute_loss_normalize(X_te, y_te, z_te, gamma, wr, nor);
                lidx = sub2ind([10, num_gamma, num_nu], i, j , k);
                ridge_results(lidx) = loss;
                
            end
        end
        
    end

    % flatten out
    simSpace = [10, num_gamma];
    numSims = prod(simSpace);
    bisect_results = zeros(numSims, 1);
    bruckner_results = zeros(numSims, 1);
    ssdp_results = zeros(numSims,1);
    socp_results = zeros(numSims,1);
    bitime = zeros(numSims, 1);
    sstime = zeros(numSims, 1);
    
    for idx=1:numSims
        [i, j] = ind2sub(simSpace, idx);
        trIdx = CVO.training(i);
        teIdx = CVO.test(i);
        X_tr = X(trIdx, :);
        y_tr = y(trIdx);
        z_tr = z(trIdx);
        X_te = X(teIdx, :);
        y_te = y(teIdx);
        z_te = z(teIdx);
        
        gamma = gamma_range(j);

        % bisection method
        fprintf('Doing bisect\n')
        [w_s0, W_s, B, Fq, time, iter] = bisect_mosek(X_tr,y_tr, z_tr, gamma, bias_upper, q_start, tol);
        bisecttime(idx) = time; 
        
        W = [W_s, w_s0; w_s0.', 1];
        w_s = w_s0(2:end);
        w_opt = rank_decompose(B, W);
        w_opt = w_opt(2:end-1);
        
        % reverse normalization loss
        bisect_results(idx)  = compute_loss_normalize(X_te, y_te, z_te, gamma, w_opt, nor);

        fprintf('Doing ssdp & socp\n');
        [w_ssdp, optval_ssdp, time_ssdp, lambda_ssdp] = singlesdp_mosek(X_tr , y_tr, z_tr, gamma);
        sstime_yalmip(idx) = time_ssdp;  
         
        [w_socp, optval_socp, time_socp, time_eig] = socp_mosek(X_tr , y_tr, z_tr, gamma);
        ssocptime_yalmip(idx) = time_socp; 
              
        ssdp_results(idx) = compute_loss_normalize(X_te, y_te, z_te, gamma, w_ssdp, nor);
        socp_results(idx) = compute_loss_normalize(X_te, y_te, z_te, gamma, w_socp, nor);

    end

    % bruckner method
    if bruckner_flag
        parfor_progress(numSims);
        parfor idx = 1:numSims
            [i, j] = ind2sub(simSpace, idx);
            trIdx = CVO.training(i);
            teIdx = CVO.test(i);
            X_tr = X(trIdx, :);
            y_tr = y(trIdx);
            z_tr = z(trIdx);
            X_te = X(teIdx, :);
            y_te = y(teIdx);
            z_te = z(teIdx);

            gamma = gamma_range(j);
            w_b = bruckner_method(X_tr, y_tr, z_tr, gamma);
            % reverse normalization loss
            bruckner_results(idx) = compute_loss_normalize(X_te, y_te, z_te, gamma, w_b, nor);
        end
    end
    
    % reshape result
    ridge_results = reshape(ridge_results, 10, num_gamma, num_nu);
    bruckner_results = reshape(bruckner_results, 10, num_gamma);
    bisect_results = reshape(bisect_results, 10, num_gamma);
    ssdp_results = reshape(ssdp_results, 10, num_gamma);
    socp_results = reshape(socp_results, 10, num_gamma);
    % get mean result
    bisect_means = mean(bisect_results, 1)';
    bruckner_means = mean(bruckner_results, 1)';
    ridge_means = mean(ridge_results, 1);
    ssdp_means = mean(ssdp_results, 1)';
    socp_means = mean(socp_results, 1)';
    ridge_best_means = max(ridge_means, [], 3)';
end






