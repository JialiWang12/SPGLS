%sdp_compare
clear all
clc
% set m, n list, where m = multiple * n
n_list = [100, 500, 1000, 2000, 4000, 6000]';
multiple_list = [0.5,1,2];
gamma_list = [1e-2, 1e-1];
% set synthetic data path
path = 'datasets/synthetic';
% set the scaling parameter of normalization 
const = 1;
% set the parameter of noise
sigma = 0.5; flag_noise = 1;

for gamma_idx = 1: length(gamma_list)
    gamma = gamma_list(gamma_idx);
    for multiple_idx = 1:length(multiple_list)
        multiple = multiple_list(multiple_idx);
        m_list = multiple * n_list;
        len_m = length(m_list);
        % initial
        bisecttime = zeros(len_m,1);ssdptime = zeros(len_m,1);socptime = zeros(len_m,1); iterlist = zeros(len_m,1);
        % flag of bisect and single SDP method
        flag_bisect = 1; flag_ssdp = 1; 
        bisect_result= zeros(len_m,1);ssdp_result = zeros(len_m,1);socp_result = zeros(len_m,1); time_eig_list = zeros(len_m,1);
        for idx = 1: len_m
                % read data
                m = m_list(idx); n = n_list(idx);
                fprintf('\ndataset size(%d,%d)\n',m,n);               
                % read dense data generated by scikit-learn function: make_regression 
                csvname = strcat(path,'/m',string(m),'n',string(n),'.csv');
                data = readtable(csvname);
                X = data{:,1:n}; y = data{:,n+1};
                % creat fake z 
                y_quantile_low = quantile(y,0.25); y_quantile_high = quantile(y,0.75);        
                z = y; z(z<y_quantile_low) = y_quantile_low; 
                % add noise
                if flag_noise
                    z = z + sigma * randn(size(z));
                end
                % normalization
                nor = const*max(abs(y));
                X = normalize(X,'range'); z = z/nor; y = y/nor;                
                
                % bisection parameter
                tol = 1e-4;q_start = 13380000 / 500.0; bias_upper = 20;
                % bisection method
                if flag_bisect 
                    try % as dimension become higher, the bisect method might have out of memory error
                        [w_s0, W_s, B, Fq, time_bisect, iter_bisect] = bisect_mosek(X, y, z, gamma, bias_upper, q_start, tol);
                        bisecttime(idx) = time_bisect; 
                        iterlist(idx)= iter_bisect;
                        % rank decompose
                        W = [W_s, w_s0; w_s0.', 1];
                        w_s = w_s0(2:end);
                        w_opt = rank_decompose(B, W);
                        w_opt = w_opt(2:end-1);
                        
                        bisect_result(idx) = compute_loss_normalize(X, y, z, gamma, w_opt, nor);
                        
                        if bisecttime(idx) > 1800
                            flag_bisect = 0;
                        end
                    catch
                        flag_bisect = 0;
                    end
                end
                
                % single SDP method
                if flag_ssdp
                    [w_ssdp,optval_ssdp,time_ssdp,lam_ssdp] = singlesdp_mosek(X, y, z, gamma);
                    ssdptime(idx) = time_ssdp; 
                    ssdp_result(idx) = compute_loss_normalize(X, y, z, gamma, w_ssdp, nor);
                    % if time of single SDP > 1800 seconds. it would not be
                    % run in higher dimension case
                    if ssdptime(idx) > 1800
                        flag_ssdp = 0;
                    end
                end

                % SOCP method 
                [w_socp,optval_socp,time_socp,time_eig] = socp_mosek(X, y, z, gamma);
                socptime(idx) = time_socp; time_eig_list(idx) = time_eig;
                socp_result(idx) = compute_loss_normalize(X, y, z, gamma, w_socp, nor);

        end
        
        % calculate the ratios
        ratio1 = bisecttime./socptime;
        ratio2 = ssdptime./socptime; 
        % sava result
        result_table = table(m_list,n_list,bisecttime,ssdptime,socptime,ratio1,ratio2,iterlist,time_eig_list);      
        table_name =strcat('./result/synthetic_result/gamma',string(gamma),'/mutiple',string(multiple),'.csv');
        writetable(result_table,table_name);
    end
    
end
















