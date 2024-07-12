addpath('../')
load('data_ori.mat')
t = data_ori.DayNumber;

figure(Position=[0, 0, 500, 250])
plot(t, data_ori.Cyclopoids, LineWidth=1, Color="red", LineStyle='-')
xlim([t(1), t(end)])
fontsize(18, "points")

figure(Position=[0, 0, 500, 250])
plot(t, data_ori.Rotifers, LineWidth=1, Color="blue", LineStyle='-')
xlim([t(1), t(end)])
fontsize(18, "points")

figure(Position=[0, 0, 500, 250])
plot(t, data_ori.Picophytoplankton, LineWidth=1, Color="#2A8947", LineStyle='-')
xlim([t(1), t(end)])
fontsize(18, "points")
