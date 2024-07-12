clc;
clear;
close all;
addpath('../../');

%% load data
load('data_ori.mat');

data = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
t = data_ori.DayNumber;

num_component = 3;

%% equidistant

t_fit = (t(1):1:t(end)).';
y_spline = zeros(length(t_fit), num_component);
y_linear = zeros(length(t_fit), num_component);
for i = 1:num_component
    y_spline(:,i) = interp1(t, data(:, i), t_fit, 'spline');
    y_linear(:,i) = interp1(t, data(:, i), t_fit, 'linear');
end

%% test GC
% parameters
num_component = 3;
sig_level = 0.05;
max_lag = 15;

% variables
cause_list = zeros(num_component,1);
F_list = zeros(num_component,1);
critic_list = zeros(num_component,1 );

%test
for i = 1:num_component
    for j = 1:1
        if i == j
            continue
        end

        [F_tmp, c_tmp] = granger_cause(y_linear(:,j), y_linear(:,i), sig_level, max_lag);
        F_list(i,j) = F_tmp;
        critic_list(i,j) = c_tmp;
        if F_tmp > c_tmp
            cause_list(i,1) = cause_list(i,1) + 1;
        end
    end
end

cause_list
