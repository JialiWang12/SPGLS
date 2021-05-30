%sdp_compare
clear all
clc


% % read data: wine_modest, wine_severe, insurance_modest,insurance_severe,
% % building_modest,building_severe,blog_modest,blog_severe,
dataname = 'wine_modest';
% add noise to fake lable z and set noise type: Guassian or Uniform or None
noise.flag = 1; noise.type = 'None'; 
% sigma represents the standard deviation of Guassin noise
noise.sigma = 0.5;
%the half width of uniform interval
noise.halfwidth = 1;
[X, y, z, const, gamma_list, gamma_time, datasize_list] = data_read(dataname,noise);


% get the results of mse or time comparison
getMSE_yes = 1; getTime_yes = 1;

% parameter for regularizer in ridge regression
nu_range = [1e-5 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100 ,1000]; 

tol = 1e-4; % tol of bisection method
if contains(dataname,'wine')
    q_start = 1000;
else
    q_start = 13380000 / 500.0;
end
bias_upper = 20;

% set bruckner flag
if contains(dataname,'blog')
    flag = 0;
else
    flag = 1;
end

if getMSE_yes
    % get mse results
    [bisect_means,bruckner_means, ridge_means, ssdp_means, socp_means] = get_mse(X, y, z, q_start, bias_upper, tol, gamma_list, nu_range, const, flag);
    % save result
    mse_table = table(gamma_list,bisect_means,bruckner_means, ridge_means, ssdp_means, socp_means);
    mse_table_name = strcat('./result/',dataname,'_',noise.type,'_mse.csv');
    writetable(mse_table,mse_table_name);
end

if getTime_yes
    % get time results
    [bruck_mr, bruck_logmr, bruck_std, bruck_logstd, bisect_mr, bisect_logmr, bisect_std, bisect_logstd, ssdp_mr,ssdp_logmr, ssdp_std, ssdp_logstd, socp_mr, socp_logmr, socp_std, socp_logstd] = get_time(X, y, z, gamma_time, q_start, tol, datasize_list, bias_upper, flag, dataname);
    % save result
    time_table = table(datasize_list,bruck_mr, bruck_logmr, bruck_std, bruck_logstd, bisect_mr, bisect_logmr, bisect_std, bisect_logstd, ssdp_mr,ssdp_logmr, ssdp_std, ssdp_logstd, socp_mr, socp_logmr, socp_std, socp_logstd);
    time_table_name = strcat('./result/',dataname,'_',noise.type,'_time.csv');
    writetable(time_table,time_table_name);
end
