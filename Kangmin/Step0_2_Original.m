clc;
clear;
close all;
addpath('../')

%% load raw data
load('rms_series_mouse-4_session-1.mat')
t = t.';
labels = ["\theta_{low}", "\theta_{high}", "\beta_{high}"];

num_component = length(y(1,:));

%% equidistant interpolation

interpolate_ratio = 1/2;

interpolated_interval = (t(2) - t(1)) * interpolate_ratio;

t_interp = (t(1):interpolated_interval:t(end)).';
y_interp = zeros(length(t_interp), num_component);
for i = 1:num_component
    y_interp(:,i) = interp1(t, y(:, i), t_interp);
end

%% plot equidistant data

data_before = y;
t_before = t;
data_after = y_interp;
t_after = t_interp;
for j = 1:num_component
    figure(i* 10 + j);
    % plot(t_before, data_before(:,j),'k',t_after, data_after(:,j),'r')
    plot(t_after,data_after(:, j), LineWidth=1);
    xlim([t(1),t(end)])
    legend(labels(j));
end

%% cut time series data
window_size_ori = 799;
overlapping_ratio = 0.9;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

y_ori = {};
start = 1;
data_tmp = y_interp;
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

num_fourier = 8;
option = ['fourier', num2str(num_fourier)];

% y_fit = y_ori;
y_fit = {};
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

%% plot interpolated data
idx = 1;

% for j = 1:7
    % idx = (j-1) * 10 + 1;
    y_ori_tmp = y_ori{idx};
    y_fit_tmp = y_fit{idx};
    for i = 1:num_component
        figure(j*10 + i)
        plot(t_fit, y_ori_tmp(:,i), Color="#0072BD", LineWidth=1);
        hold on
        plot(t_fit, y_fit_tmp(:,i), Color='r', LineWidth=2);
        % legend(labels(i))
    end
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
% y = data_ori;
y_total = y_fit;
time_interval = 1/window_size_ori;

save('data_cut', 't', 'y', 'y_total', 'time_interval', 'num_data', 'num_component', 'labels')
