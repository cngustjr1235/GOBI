clc;
clear;
close all;

%% load raw data
load('raw_data_repressilator.mat')
t_ori = t;
y_ori = y;
labels = ["TetR", "Lacl", "\lambda cl"];

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
	y_fit(:,i) = interp1(t_ori, y_ori(:, i), t_fit, 'spline');
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

tetr = y_fit(:, 1);
lacl = y_fit(:, 2);
lcl = y_fit(:, 3);
tetrdot = gradient(tetr, time_interval);
lacldot = gradient(lacl, time_interval);
lcldot = gradient(lcl, time_interval);

fixed_target = 4:6;

y_ext = [tetr, lacl, lcl, tetrdot, lacldot, lcldot];

[y_total, t_total, class_pair_list, class_count_list] = cut_horizontal(y_ext, t_fit, {[], [], [], [], [], []});

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
