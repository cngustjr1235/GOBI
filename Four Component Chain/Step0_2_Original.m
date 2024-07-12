clc;
clear;
close all;
addpath('../')

%% load raw data
load('data.mat')
labels = ["y1", "y2", "y3", "y4"];

num_component = length(labels);
length_timeseries = length(y(:,1));

%% equidistant interpolation
t_fit = (t(1):0.5:t(end)).';
data_fit = zeros(length(t_fit), num_component);
for i = 1:num_component
    data_fit(:,i) = interp1(t, y(:, i), t_fit, 'spline');
end

%% cut time series data
window_size_ori = 30;
overlapping_ratio = 0.1;

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

%% interpolate data


% num_fourier = 8;
% option = ['fourier', num2str(num_fourier)];
% 
% y_fit = {};
y_fit = y_ori;
t_fit = linspace(0, 1, window_size_ori).';
% 
% for i = 1:num_data
%     y_fit_tmp = y_ori{i};
%     for j = 1:num_component
%         fitting = fit(t_fit, y_fit_tmp(:, j), option);
%         w = fitting.w;
%         fouriers = ones(1, length(t_fit));
%         for k = 1:num_fourier
%             fouriers = [fouriers; cos(k*w*t_fit.')];
%             fouriers = [fouriers; sin(k*w*t_fit.')];
%         end
%         coeffs = coeffvalues(fitting);
%         y_fit_tmp(:,j) = coeffs(1:end-1) * fouriers;
%     end
%     y_fit{end+1} = y_fit_tmp;
% end

% %% plot interpolated data
% idx = 1;
% 
% for j = 1:7
%     idx = (j-1) * 10 + 1;
%     y_ori_tmp = y_ori{idx};
%     y_fit_tmp = y_fit{idx};
%     for i = 1:num_component
%         figure(j*10 + i)
%         plot(t_fit, y_ori_tmp(:,i), 'k', t_fit, y_fit_tmp(:,i), 'r-o');
%         legend(labels(i))
%     end
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
y_total = y_fit;
time_interval = 1/window_size_ori;

save('data_cut', 't', 'y', 'y_total', 'time_interval', 'num_data', 'num_component', 'labels')
