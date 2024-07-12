clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
t = data_ori.DayNumber;
data_ori = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
data_ori = data_ori(1:end,:);
t = t(1:end, :);
labels = ["Cyclopoids", "Protozoa", "Rotifers"];

num_component = length(labels);

data_ori(isnan(data_ori)) = 0;

%% equidistant interpolation by smooth splining
t_fit = (t(1):1:t(end)).';
data_fit = zeros(length(t_fit), num_component);
for i = 1:num_component
    data_fit(:,i) = csaps(t, data_ori(:, i), 0.01, t_fit);
end

%% plot equidistant data

data_before = data_ori;
t_before = t;
data_after = data_fit;
t_after = t_fit;
for j = 1:num_component
    figure(i* 10 + j);
    plot(t_before, data_before(:,j),'k',t_after, data_after(:,j),'r')
    xlim([t(1),t(end)])
    legend(labels(j));
end

%% cut time series data
window_size_ori = 60;
overlapping_ratio = 0.9;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

y_max = zeros(num_component, 1);
y_min = zeros(num_component, 1);
for i = 1:num_component
    y_max(i) = max(data_fit(:,i));
    y_min(i) = min(data_fit(:,i));
end

data_precut = {data_fit(216:792, :), data_fit(1027:end, :)};
y_ori = {};

for i = 1:length(data_precut)
    start = 1;
    data_tmp = data_precut{i};
    length_timeseries = length(data_tmp(:,1));
    while(1)
        if start + window_size_ori > length_timeseries
            start = length_timeseries - window_size_ori + 1;
        end
        y_tmp = data_tmp(start:start+window_size_ori-1, :);
        for j = 1:num_component
            % normalize data
            y_tmp(:,j) = y_tmp(:,j) - y_min(j);
            if max(y_tmp(:,j)) ~= 0
                y_tmp(:,j) = y_tmp(:,j) / y_max(j);
            end
        end
        y_ori{end+1} = y_tmp;
        if start + window_size_ori > length_timeseries
            break;
        end
        start = start + window_move_ori;
    end
end

num_data = length(y_ori);

% %% plot cutted data
% for i = 1:1
%     data_tmp = y_ori{i};
%     t_tmp = linspace(0,1,window_size_ori);
%     figure(i);
%     for j = 1:num_component
%         plot(t_tmp, data_tmp(:,j), LineWidth=1);
%         hold on
%     end
% end

%% interpolate data

y_fit = y_ori;
t_fit = linspace(0, 1, window_size_ori).';

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
