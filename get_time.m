function [bruck_mr, bruck_logmr, bruck_std, bruck_logstd, bisect_mr, bisect_logmr, bisect_std, bisect_logstd, ssdp_mr,ssdp_logmr, ssdp_std, ssdp_logstd, socp_mr, socp_logmr, socp_std, socp_logstd] = get_time(X, y, z, gamma, q_start, tol, size_range, bias_upper,varargin)
    num_sizes = length(size_range);
    
    bruck_mr = zeros(num_sizes, 1);
    bruck_logmr = zeros(num_sizes, 1);
    bruck_std = zeros(num_sizes, 1);
    bruck_logstd = zeros(num_sizes, 1);
    
    bisect_mr = zeros(num_sizes, 1);
    bisect_logmr = zeros(num_sizes, 1);
    bisect_std = zeros(num_sizes, 1);
    bisect_logstd = zeros(num_sizes, 1);
    
    ssdp_mr = zeros(num_sizes, 1);
    ssdp_logmr = zeros(num_sizes, 1);
    ssdp_std = zeros(num_sizes, 1);
    ssdp_logstd = zeros(num_sizes, 1);
    
    socp_mr = zeros(num_sizes, 1);
    socp_logmr = zeros(num_sizes, 1);
    socp_std = zeros(num_sizes, 1);
    socp_logstd = zeros(num_sizes, 1);  
    
    if nargin == 9
         bruckner_flag = varargin{1};
    end
    
    if nargin == 10
         bruckner_flag = varargin{1};
         dataname = varargin{2};
    end

    fold = 10;

    for i=1:num_sizes   
        m = size_range(i);
        fprintf('Starting a dataset size m = %d\n',m);
        bruck_m = zeros(1, fold);
        bisect_m = zeros(1, fold);
        ssdp_m = zeros(1, fold);
        socp_m = zeros(1, fold);
        
        bruck_fval_m = zeros(1, fold);
        bisect_fval_m = zeros(1, fold);
        ssdp_fval_m = zeros(1, fold);
        socp_fval_m = zeros(1, fold);
        
        bruck_iter_m = zeros(1, fold);
        for j=1:fold
            perm = randperm(size(X, 1));
            X_shuffle = X(perm, :);
            y_shuffle = y(perm, :);
            z_shuffle = z(perm, :);
            
            X_tr = X_shuffle(1:m , :);
            y_tr = y_shuffle(1:m);
            z_tr = z_shuffle(1:m);
            
            [w_s0, W_s, B, Fq, time, iter] = bisect_mosek(X_tr,y_tr, z_tr, gamma, bias_upper, q_start, tol);
            bisect_m(1,j) = time;
            
            
            W = [W_s, w_s0; w_s0.', 1];
            w_s = w_s0(2:end);
            w_opt = rank_decompose(B, W);
            w_opt = w_opt(2:end-1);
            bisect_fval_m(1,j) = fval(X_tr,y_tr,z_tr,gamma,w_opt);
            
            tic;
            if bruckner_flag
                [w_br, iter_br] = bruckner_method(X_tr, y_tr, z_tr, gamma);
                bruck_m(1,j) = toc;
                bruck_fval_m(1,j) = fval(X_tr,y_tr,z_tr,gamma,w_br);
                bruck_iter_m(1,j) = iter_br;
            end
            
            [w_socp, ~, time_socp,~] = socp_mosek(X_tr,y_tr,z_tr,gamma);
            socp_m(1,j) = time_socp;
            socp_fval_m(1,j) = fval(X_tr,y_tr,z_tr,gamma,w_socp);
            
            [w_ssdp, ~, time_ssdp,~] = singlesdp_mosek(X_tr,y_tr, z_tr, gamma);
            ssdp_m(1,j) = time_ssdp;
            ssdp_fval_m(1,j) = fval(X_tr,y_tr,z_tr,gamma,w_ssdp);

           

        end
        
        bruck_mr(i, 1) = mean(bruck_m);
        bruck_logmr(i, 1) = mean(log10(bruck_m));
        bruck_std(i, 1) = std(bruck_m);
        bruck_logstd(i, 1) = std(log10(bruck_m));
        
        bisect_mr(i, 1) = mean(bisect_m);
        bisect_logmr(i, 1) = mean(log10(bisect_m));
        bisect_std(i, 1) = std(bisect_m);
        bisect_logstd(i, 1) = std(log10(bisect_m));
        
        ssdp_mr(i,1) = mean(ssdp_m);
        ssdp_logmr(i, 1) = mean(log10(ssdp_m));
        ssdp_std(i,1) = std(ssdp_m);
        ssdp_logstd(i,1) = std(log10(ssdp_m));
        
        socp_mr(i,1) = mean(socp_m);
        socp_logmr(i, 1) = mean(log10(socp_m));
        socp_std(i,1) = std(socp_m);
        socp_logstd(i,1) = std(log10(socp_m));
        
    end
 
end