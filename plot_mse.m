function plot_mse(csvname)
table = readtable(csvname);
varname_list = table.Properties.VariableNames;
for i = 1: length(varname_list)
    eval( strcat(string(varname_list(i)),'=table{:,i};'))
end
fig = figure('Name',csvname);
figure(1);
hold on
if sum(bruckner_means) ~= 0 % it means we do not gather bruckner result
    plot(gamma_list, bruckner_means, 'LineStyle',  '-.',  'color', 'blue','linewidth',2);
end
plot(gamma_list, ridge_means, 'LineStyle', '-', 'color', 'red','linewidth',2);
plot(gamma_list, bisect_means,'LineStyle', '-', 'color', 'g','linewidth',2);
plot(gamma_list, ssdp_means, 'LineStyle', '-.', 'color', 'm','linewidth',3);
plot(gamma_list, socp_means, 'LineStyle',  ':',  'color', 'k','linewidth',2);

if sum(bruckner_means) == 0 % it means we do not gather bruckner result
    legend('Ridge Regression', 'Bisection SDP','Single SDP', 'SOCP')
else
    legend('Br$\mathrm{\ddot{u}}$ckner and Scheffer (2011)','Ridge Regression','Bisection SDP', 'Single SDP', 'SOCP','Interpreter','latex')
end

xlabel('$\gamma$ ','Interpreter','latex');
ylabel('MSE (modest)');
legend('Location','east','FontSize',10);

hold off



