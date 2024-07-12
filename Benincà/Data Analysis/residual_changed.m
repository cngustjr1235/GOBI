% To calculate residual of initial interpolation
% Cut, Normalize, interpolate with Fourier8

clc;
clear;
close all;

%% import data
load('data.mat')
% t = data.DayNumber;
y = horzcat(data.Rotifers, data.Calanoids, data.Picophytoplankton, data.Nanophytoplankton);

%% parameters
num_component = length(y(1,:));
window_size_ori = 140;
overlapping_ratio = 0.1;
time_interval = 1/10;

t = linspace(0, window_size_ori-1,window_size_ori).';

%% process parameters
window_size = window_size_ori;
window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));
window_move = window_move_ori;

%% cut data
disp('cut_data...')
y_total = {};
start = 1;
length_timeseries = length(y(:,1));
while(1)
    if start + window_size > length_timeseries
        break
    end
    y_tmp = y(start:start + window_size - 1, :);
    for i = 1:num_component
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i));
    end
    y_total{end+1} = y_tmp;
    start = start + window_move;
end
length_data = length(y_total);

%% plot cutted data
for j = 1:length_data
    figure(j);
    y_tmp = y_total{j};
    for i = 1:num_component
        plot(t, y_tmp(:,i));
        hold on;
    end
end

%% interpolate
% interpolation method
% method = 1: linear interpolation
% method = 2: spline interpolation
% method = 3: fourier interpolation, 
%             In this case, users have to choose the order of fourier method (1~8)

method = 3;
num_fourier = 8;

y_fit_total = {};
t = linspace(0, window_size_ori-1,window_size_ori).';
t_fit = linspace(0, window_size_ori-1, (window_size_ori-1)/time_interval+1).';

for i = 1:length_data
    y_fit = zeros(length(t_fit), num_component);
    y_tmp = cell2mat(y_total(i));
    for j = 1:num_component
        y_fit(:,j) = interp1(t,y_tmp(:,j),t_fit,'linear');
    end
    option = ['fourier', num2str(num_fourier)];
    for j = 1:num_component
        fitting = fit(t_fit, y_fit(:,j), option);
        w = fitting.w;
        fouriers = ones(1,length(t_fit));
        for k = 1:num_fourier
            fouriers = [fouriers; cos(k*w*t_fit.')];
            fouriers = [fouriers; sin(k*w*t_fit.')];
        end
        coeffs = coeffvalues(fitting);
        y_fit(:,j) = coeffs(1:end-1) * fouriers;
    end
    y_fit(y_fit < 0) = 0;
    y_fit_total{end+1} = y_fit;
end

%% plot interpolated data
for j = 1:length_data
    figure(j);
    y_tmp = y_total{j};
    y_fit_tmp = y_fit_total{j};
    for i = 1:num_component
        plot(t, y_tmp(:,i), t_fit, y_fit_tmp(:,i));
        hold on
    end
end

%% calculate residual
residual_total = zeros(length_data, num_component);
for i = 1:length_data
    y_fit = cell2mat(y_fit_total(i));
    y_fit = y_fit(1:1/time_interval:end, :);
    y_noise = cell2mat(y_total(i));
    for j = 1:num_component
        residual_tmp = mean((y_fit(:,j)- y_noise(:,j)).^2);
        residual_total(i,j) = residual_tmp;
    end
end
mean(mean(residual_total))
