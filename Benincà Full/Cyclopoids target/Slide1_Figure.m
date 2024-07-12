addpath('../')
load('data_ori.mat')
t = data_ori.DayNumber;
y = data_ori.Rotifers;
% y(y>60) = nan;

figure(Position=[0, 0, 1200, 400])

plot(t, y, LineWidth=1, Color="blue", LineStyle='-')

xlim([t(1), t(end)])
xlabel("Time(Days)")
ylabel("Rotifers")
fontsize(18, "point")

ylim([0,5])