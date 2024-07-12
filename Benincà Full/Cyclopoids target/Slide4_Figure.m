clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
y = data_ori.Protozoa;
y(isnan(y)) = 0;
t = data_ori.DayNumber;
% y = nthroot(y, 4);

%% equidistant interpolation

t_fit = (t(1):1:t(end)).';
y_fit = interp1(t, y, t_fit, 'linear');

t_cut1 = t_fit(361:480);
y_cut1 = y_fit(361:480);

% t_cut2 = t_fit(2539:2658);
% y_cut2 = y_fit(2539:2658);

y_cut1 = y_cut1 - min(y_cut1);
if max(y_cut1) ~= 0
    y_cut1 = y_cut1 / max(y_cut1);
end

% y_cut2 = y_cut2 - min(y_cut2);
% y_cut2 = y_cut2 / max(y_cut2);

y_fit1 = y_cut1;
% y_fit2 = y_cut2;
t_fit = t_cut1 / t_cut1(end);

figure(Position=[0, 0, 500, 400])
plot(t_fit, y_fit1, LineWidth=2, Color="#2A8947", LineStyle="-")
xlim([t_fit(1), t_fit(end)])
fontsize(20, "points")
yticks([0, 1])
xticks([0, 1])

% figure(Position=[0, 0, 500, 400])
% plot(t_fit, y_fit2, LineWidth=2, Color="blue", LineStyle="-")
% xlim([t_fit(1), t_fit(end)])
% fontsize(20, "points")
% yticks([0, 1])
% xticks([0, 1])

num_fourier = 8;
option = ['fourier', num2str(num_fourier)];

fitting = fit(t_fit, y_fit1, option);
w = fitting.w;
fouriers = ones(1, length(t_fit));
for k = 1:num_fourier
    fouriers = [fouriers; cos(k*w*t_fit.')];
    fouriers = [fouriers; sin(k*w*t_fit.')];
end
coeffs = coeffvalues(fitting);
y_fit1 = coeffs(1:end-1) * fouriers;
% 
% fitting = fit(t_fit, y_fit2, option);
% w = fitting.w;
% fouriers = ones(1, length(t_fit));
% for k = 1:num_fourier
%     fouriers = [fouriers; cos(k*w*t_fit.')];
%     fouriers = [fouriers; sin(k*w*t_fit.')];
% end
% coeffs = coeffvalues(fitting);
% y_fit2 = coeffs(1:end-1) * fouriers;


figure(Position=[0, 0, 500, 400])
plot(t_fit, y_cut1, LineWidth=2, Color="#2A8947", LineStyle="-")
hold on
plot(t_fit, y_fit1, LineWidth=3, Color="#D95319", LineStyle="--")
xlim([t_fit(1), t_fit(end)])
ylim([0,1])
fontsize(20, "points")
yticks([0, 1])
xticks([0, 1])

% figure(Position=[0, 0, 500, 400])
% plot(t_fit, y_cut2, LineWidth=1, Color="black", LineStyle="-")
% hold on
% plot(t_fit, y_fit2, LineWidth=2, Color="blue", LineStyle="-")
% xlim([t_fit(1), t_fit(end)])
% ylim([0,1])
% fontsize(20, "points")
% yticks([0, 1])
% xticks([0, 1])
