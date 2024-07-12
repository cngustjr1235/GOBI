clc;
clear;
close all;
addpath('../../');

%% load data
load('data_ori.mat');

data = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
t = data_ori.DayNumber;

num_component = 3;

% load('data_cardio.mat');
% y = y(1:720, :);
% tau = mdDelay(y, 'maxLag', 25, 'plottype', 'all');

%% equidistant

t_fit = (t(1):1:t(end)).';
y_spline = zeros(length(t_fit), num_component);
y_linear = zeros(length(t_fit), num_component);
for i = 1:num_component
    y_spline(:,i) = interp1(t, data(:, i), t_fit, 'spline');
    y_linear(:,i) = interp1(t, data(:, i), t_fit, 'linear');
end

% tau_spline = mdDelay(y_spline, 'maxLag', 30, 'plottype', 'all');
tau_linear = mdDelay(y_linear, 'maxLag', 30, 'plottype', 'all');