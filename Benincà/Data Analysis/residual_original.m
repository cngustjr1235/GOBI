% To calculate residual of initial interpolation
% Compress 1/4, interpolate with Fourier8, and then cut
% (while figures do cut, normalize, interpolate)

clc;
clear;
close all;

%% load raw data
load('data.mat')
t = data.DayNumber;
y_raw = horzcat(data.Rotifers, data.Calanoids, data.Picophytoplankton, data.Nanophytoplankton);
y = nthroot(y_raw, 4);

%% parameters
num_component = length(y(1,:));
window_size_ori = 90;
overlapping_ratio = 0.1;
time_interval = 1/10;


%% process parameters
window_size = window_size_ori / time_interval;
window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));
window_move = window_move_ori / time_interval;

%% interpolate
% interpolation method
% method = 1: linear interpolation
% method = 2: spline interpolation
% method = 3: fourier interpolation, 
%             In this case, users have to choose the order of fourier method (1~8)

method = 3;
num_fourier = 8;

t_fit = linspace(0,t(end),length(t)/time_interval+1).';
y_fit = zeros(length(t_fit),num_component);

for i = 1:num_component
    y_fit(:,i) = interp1(t, y(:,i),t_fit,'linear');
end
option = ['fourier',num2str(num_fourier)];
for i = 1:num_component
    fitting = fit(t_fit,y_fit(:,i),option);
    w = fitting.w;
    fouriers = ones(1,length(t_fit));
    for j = 1:num_fourier
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_fit(:,i) = coeffs(1:end-1) * fouriers;
end
y_fit(y_fit < 0) = 0;

%% plot interpolated data
for i = 1:num_component
    figure(i);
    plot(t, y(:,i), t_fit, y_fit(:,i));
end

%% cut the data
y_total = {};
y_total_ori = {};
start = 1; % it can be changed if users do not want to use the beginning of time series
start_ori = 1;
length_timeseries = length(y_fit(:,1)); % it can be changed if users do not want to use the end of time series
while(1)
    if start + window_size-1 > length_timeseries
        break
    end
    
    % cut
    y_tmp = y_fit(start:start + window_size - 1,:);
    y_tmp_ori = y(start_ori:start_ori + window_size_ori - 1,:);
    % t = t_fit(1:window_size) / t_fit(window_size);
   
    % normalize
    for i = 1:num_component
        min_ori = min(y_tmp_ori(:, i));
        max_ori = max(y_tmp_ori(:, i));
        y_tmp(:,i) = y_tmp(:,i) - min_ori;
        y_tmp(:,i) = y_tmp(:,i) / max_ori;
        y_tmp_ori(:,i) = y_tmp_ori(:,i) - min_ori;
        y_tmp_ori(:,i) = y_tmp_ori(:,i) / max_ori;
    end
    y_total{end+1} = y_tmp;
    y_total_ori{end+1} = y_tmp_ori;
    start = start + window_move;
    start_ori = start_ori + window_move_ori;
end
num_data = length(y_total);

%% plot cutted data
t = linspace(0, window_size_ori-1,window_size_ori).';
t_fit = linspace(0, window_size_ori-1, window_size_ori/time_interval).';

for j = 1:num_data
    figure(j);
    y_tmp = y_total{j};
    y_ori_tmp = y_total_ori{j};
    for i = 1:num_component
        plot(t_fit, y_tmp(:,i), t, y_ori_tmp(:,i));
        hold on
    end
end

%% calculate residual
residual_total = zeros(num_data, num_component);
for i = 1:num_data
    y_fit_tmp = cell2mat(y_total(i));
    y_fit = y_fit_tmp(1:(1/time_interval):end, :);
    y_noise = cell2mat(y_total_ori(i));
    for j = 1:num_component
        residual_tmp = mean((y_fit(:, j) - y_noise(:, j)).^2);
        residual_total(i,j) = residual_tmp;
    end
end
mean(mean(residual_total))