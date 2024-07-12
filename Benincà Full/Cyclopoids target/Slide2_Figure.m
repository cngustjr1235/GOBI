addpath('../')
load('data_ori.mat')
t = data_ori.DayNumber;

t = t(539:581);
data_ori = data_ori(539:581,:);
data_ori.Protozoa(isnan(data_ori.Protozoa)) = 0;

figure(Position=[0, 0, 1000, 200])
plot(t, data_ori.Cyclopoids, LineWidth=1, Color="red", LineStyle='-')
xlim([t(1), t(end)])
fontsize(18, "points")

hold on
plot(t, data_ori.Protozoa, LineWidth=1, Color="blue", LineStyle='-')

legend(["Cyclopoids", "Protozoa"]);

% 
% figure(Position=[0, 0, 500, 250])
% plot(t, data_ori.Rotifers, LineWidth=1, Color="blue", LineStyle='-')
% xlim([t(1), t(end)])
% fontsize(18, "points")
% 
% figure(Position=[0, 0, 500, 250])
% plot(t, data_ori.Protozoa, LineWidth=1, Color="#2A8947", LineStyle='-')
% xlim([t(1), t(end)])
% fontsize(18, "points")
