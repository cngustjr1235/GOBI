clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_cardio.mat')
t = linspace(0,1031,1032).';
data_ori = [cardio, no2, so2, o3, rspar];
labels = ["cardio", "no2", "so2", "o3", "rspar"];
num_component = 5;

%% suppress data by 1/4
data_suppress = data_ori;
data_suppress(isnan(data_suppress)) = 0;

%% equidistant interpolation
% t_fit = (t(1):1:t(end)).';
% data_fit = zeros(length(t_fit), num_component);
% for i = 1:num_component
%     data_fit(:,i) = interp1(t, data_suppress(:, i), t_fit, 'linear');
% end
% 
% %% plot equidistant data
% 
% data_before = data_suppress;
% t_before = t;
% data_after = data_fit;
% t_after = t_fit;
% for j = 1:num_component
%     figure(i* 10 + j);
%     plot(t_before, data_before(:,j), 'k',t_after, data_after(:,j), 'r')
%     legend(labels(j));
% end

%% moving average
win_avg = 7;

data_fit = data_ori;
for i = 1:num_component
    data_tmp = data_fit(:, i);
    data_fit(:,i) = movmean(data_tmp, win_avg);
end

%% interpolation
time_interval = 1/2;

t_smooth = linspace(0, 1031, 1031/time_interval+1).';
y_smooth = zeros(length(t_smooth),5);
y_tmp = data_fit;
for j = 1:5
    y_fit_tmp = interp1(t,y_tmp(:,j),t_smooth);
    y_smooth(:,j) = y_fit_tmp;
end

%% cut time series data
window_size_ori = 365 * 2;
overlapping_ratio = (365-30)/365;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

y_ori = {};
start = 1;
data_tmp = data_fit;
length_timeseries = length(data_tmp(:,1));
while(1)
    if start + window_size_ori > length_timeseries
        start = length_timeseries - window_size_ori + 1;
    end
    y_tmp = data_tmp(start:start+window_size_ori-1, :);
    for j = 1:num_component
        % normalize data
        y_tmp(:,j) = y_tmp(:,j) - min(y_tmp(:,j));
        if max(y_tmp(:,j)) ~= 0
            y_tmp(:,j) = y_tmp(:,j) / max(y_tmp(:,j));
        end
    end
    y_ori{end+1} = y_tmp;
    if start + window_size_ori > length_timeseries
        break;
    end
    start = start + window_move_ori;
end

num_data = length(y_ori);

%% plot cutted data
for i = 1:3
    data_tmp = y_ori{i};
    t_tmp = linspace(0,1,window_size_ori);
    figure(i);
    for j = 1:num_component
        plot(t_tmp, data_tmp(:,j), LineWidth=1);
        hold on
    end
end

%% interpolate data

num_fourier = 4;
option = ['fourier', num2str(num_fourier)];

y_fit = {};
% y_fit = y_ori;
t_fit = linspace(0, 1, window_size_ori).';

for i = 1:num_data
    y_fit_tmp = y_ori{i};
    for j = 1:num_component
        fitting = fit(t_fit, y_fit_tmp(:, j), option);
        w = fitting.w;
        fouriers = ones(1, length(t_fit));
        for k = 1:num_fourier
            fouriers = [fouriers; cos(k*w*t_fit.')];
            fouriers = [fouriers; sin(k*w*t_fit.')];
        end
        coeffs = coeffvalues(fitting);
        y_fit_tmp(:,j) = coeffs(1:end-1) * fouriers;
    end
    y_fit{end+1} = y_fit_tmp;
end

% %% plot interpolated data
% idx = 100;
% 
% y_ori_tmp = y_ori{idx};
% y_fit_tmp = y_fit{idx};
% for i = 1:num_component
%     figure(i)
%     plot(t_fit, y_ori_tmp(:,i), 'k', t_fit, y_fit_tmp(:,i), 'r-o');
%     legend(labels(i))
% end

%% calculate residual
residual_total = zeros(num_data, num_component);
for i = 1:num_data
    y_fit_tmp = y_fit{i};
    y_ori_tmp = y_ori{i};
    for j = 1:num_component
        residual_tmp = mean((y_fit_tmp(:,j)-y_ori_tmp(:,j)).^2);
        residual_total(i, j) = residual_tmp;
    end
end
mean(mean(residual_total))

%% save data
t = t_fit;
y = data_ori;
y_total = y_fit;
time_interval = 1/window_size_ori;

save('data_cut', 't', 'y', 'y_total', 'time_interval', 'num_data', 'num_component', 'labels')
