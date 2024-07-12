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
y_fit = interp1(t, y, t_fit, ['linear']);

t_cut = t_fit(361:480);
y_cut = y_fit(361:480);

figure(Position=[0, 0, 500, 400])
plot(t_cut, y_cut, LineWidth=2, Color="#2A8947", LineStyle="-")
xlim([t_cut(1), t_cut(120)])
fontsize(20, "points")
% 
% figure(Position=[0,0,500,400])
% plot(t_fit(2539:2658), y_fit(2539:2658), LineWidth=2, Color="blue", LineStyle="-")
% xlim([t_fit(2539), t_fit(2658)])
% fontsize(20, "points")
