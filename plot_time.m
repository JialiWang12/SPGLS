function plot_time(csvname)
table = readtable(csvname);
varname_list = table.Properties.VariableNames;
for i = 1: length(varname_list)
    eval( strcat(string(varname_list(i)),'=table{:,i};'))
end
bruck_conf = zeros(size(datasize_list, 1), 1);
bisect_conf = zeros(size(datasize_list, 1), 1);
ssdp_conf = zeros(size(datasize_list, 1), 1);
socp_conf = zeros(size(datasize_list, 1), 1);

N = 10;
SEM_bruck = bruck_logstd/ sqrt(N);
SEM_bisect = bisect_logstd / sqrt(N);
SEM_ssdp = ssdp_logstd/ sqrt(N);
SEM_socp = socp_logstd/ sqrt(N);
CI95 = tinv(0.975, N-1);

for i=1:size(bruck_mr, 1)
   bruck_conf(i, 1) = SEM_bruck(i,1) * CI95;
   bisect_conf(i, 1) = SEM_bisect(i, 1) * CI95;
   ssdp_conf(i, 1) = SEM_ssdp(i, 1) * CI95;
   socp_conf(i, 1) = SEM_socp(i, 1) * CI95;
end

figure(1)
hold on
xlabel('m');
ylabel('Running time (seconds)');
if sum(bruck_mr) ~= 0  % it means we gather bruckner result
    errorbar(datasize_list, bruck_logmr, bruck_conf, 'LineStyle',  '-.','LineWidth',2,'color',[0,0.4470,0.7410]);
end
errorbar(datasize_list, bisect_logmr, bisect_conf, 'LineStyle', '-.','LineWidth',2,'color',[0.8500,0.3250,0.0980]);
errorbar(datasize_list, ssdp_logmr, ssdp_conf, 'LineStyle',  '-.' ,'LineWidth',2,'color',[0.9290,0.6940,0.1250]);
errorbar(datasize_list, socp_logmr, socp_conf, 'LineStyle',  '-.' ,'LineWidth',2,'color',[0.4940,0.1840,0.5560]);
% set(gca,'YScale','log');
set(0,'defaultTextInterpreter','latex');
if sum(bruck_mr) ~= 0  % it means we gather bruckner result
    legend('Br$\mathrm{\ddot{u}}$ckner and Scheffer (2011)','Bisection SDP','Single SDP','SOCP','Location','best','Interpreter','latex')  
    max_y = (int32(max(bruck_mr)/5) + 1) * 5;  
else
    legend('Bisection SDP','Single SDP','SOCP','Location','east')
    max_y = (int32(max(bisect_mr)/5) + 1) * 5;    
end
xlim([min(datasize_list),max(datasize_list)]);
yTick = unique(int32(get(gca,'ytick')));
set(gca,'YTick',[min(yTick):1:max(yTick)])
yTickLabels = cellstr(num2str((yTick(:)),'10^{%d}'));
yticklabels(yTickLabels);

hold off
