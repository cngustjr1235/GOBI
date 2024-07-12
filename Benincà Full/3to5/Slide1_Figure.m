addpath('../')
load('data_nut_ori.mat')
t = data_nut_ori.DayNumber;

figure(Position=[0, 0, 500, 250])
colororder({'k','k'})

yyaxis left
plot(t, data_nut_ori.TotalDissolvedInorganicNitrogen, LineWidth=1, Color="red", LineStyle='-')
hold on
% plot(t, data_ori.Picophytoplankton, LineWidth=1, Color="black", LineStyle='-')


yyaxis right
plot(t, data_nut_ori.SolubleReactivePhosphorus, LineWidth=1, Color="black", LineStyle='-')

xlim([t(1), t(end)])
xticks([])

yticks([])
yyaxis left
yticks([])