clc;
clear;
close all;

%% load raw data
load('data_ori.mat')
t_ori = data_ori.DayNumber;
y_ori = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
labels = ["Cyclopoids", "Protozoa", "Rotifers"];

num_component = length(y_ori(1,:));
y_ori(isnan(y_ori)) = 0;

%% parameters - cut the time series data
window_size_ori = t_ori(end);     % For oscillatory data, 1 period is recommended
overlapping_ratio = 0.1;  % overlapping ratio of moving window technique

% choose sampling rate for interpolation. 
% window_size_ori/time_interval becomes number of time points per cutted data
% window_size_ori/time_interval = 100 is recommended
% window_size_ori/time_interval is high (low) make the inference accurate (less accurate) and slow (fast)
time_interval = 1;       

%% process parameters
window_size = window_size_ori / time_interval;
window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));
window_move = window_move_ori / time_interval;


%% plot raw data
for i = 1:num_component
    figure(i)
    plot(t_ori, y_ori(:,i))
    xlim([0,t_ori(end)])
    xticks([])
    ylim([0,max(y_ori(:,i))])
    yticks([])
end

%% interpolate data
spline_param = 0.01;

t_fit = (t_ori(1):time_interval:t_ori(end)).';
y_fit = zeros(length(t_fit),num_component);

for i = 1:num_component
	y_fit(:,i) = csaps(t_ori, y_ori(:, i), spline_param, t_fit);
end

%% plot interpolated data
for i = 1:num_component
    figure(i)
    plot(t_fit, y_fit(:,i))
    xlim([0,t_fit(end)])
    xticks([])
    ylim([0,max(y_fit(:,i))])
    yticks([])
end

%% cut and interpolate the data
y_total = {};
start = 1; % it can be changed if users do not want to use the beginning of time series
length_timeseries = length(y_fit(:,1)); % it can be changed if users do not want to use the end of time series
while(1)
    if start + window_size-1 > length_timeseries
        break
    end
    
    % cut
    y_tmp = y_fit(start:start + window_size - 1,:);
    t = t_fit(1:window_size) / t_fit(window_size);
   
    % normalize
    for i = 1:num_component
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i));
    end
    y_total{end+1} = y_tmp; 
    start = start + window_move;
end
num_data = length(y_total);

%% plot cutted data
figure(3)
for j = 1:num_data
    y_tmp = cell2mat(y_total(j));
    for i = 1:num_component
        plot(t, y_tmp(:,i))
        hold on
    end
end
xlim([0,t(end)])

xlim([0,t(end)])
xticks([])
ylim([0,1.02])
yticks([])
%% Save data
time_interval = 1/window_size;
save('data_cut','t_fit','y_fit','y_total','time_interval','num_data','num_component')
