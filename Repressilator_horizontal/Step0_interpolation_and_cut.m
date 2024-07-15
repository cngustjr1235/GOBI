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

% choose sampling rate for interpolation. 
time_interval = 1;       

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

%% cut the data

cyc = y_fit(:, 1);
pro = y_fit(:, 2);
rot = y_fit(:, 3);
cycdot = gradient(cyc, time_interval);

fixed_target = 4;

y_ext = [cyc, pro, rot, cycdot];

[y_total, t_total, class_pair_list, class_count_list] = cut_horizontal(y_ext, t_fit, {[], [], [], []});

%% normalize
num_data = length(y_total);

for idx_data = 1:num_data
    y_data = y_total{idx_data};
    for idx_component = 1:num_component
	y_comp = y_data(:, idx_component);
	y_comp = y_comp - min(y_comp);
	y_comp = y_comp / max(y_comp);
	y_data(:, idx_component) = y_comp;
    end
    y_total{idx_data} = y_data;
end

%% Save data

save('data_cut','t_fit','y_fit','t_total','y_total','fixed_target','time_interval','num_data','num_component')
