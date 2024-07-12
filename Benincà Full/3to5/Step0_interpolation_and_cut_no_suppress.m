clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
data_ori = removevars(data_ori, {'Bacteria', 'FilamentousDiatoms','Harpacticoids', 'Ostracods', 'CalanoidCopepods', 'Nanophytoplankton', 'Protozoa' });
labels = string(data_ori.Properties.VariableNames);
labels = labels(2:end);
num_component = length(labels);

%% filtering long zero & nan(Cyclopoids, Rotifers, Protozoa) values
% Filtering for cyclopoids
% data_filtered = {data_ori(78:113, :),data_ori(126:165, :), data_ori(489:552, :), data_ori(591:666, :), data_ori(692:756, :)};

% Filtering for Calanoids

% data_filtered = {data_ori(1:206, :), data_ori(249:391, :), data_ori(488:552, :), data_ori(658:695, :)};
% data_filtered = {data_ori(1:255, :), data_ori(340:552, :), data_ori(591:end, :)};
data_filtered = data_ori;

%% suppress data by 1/4
data_suppress = {data_filtered{:,2:end}};
t_suppress = {data_filtered.DayNumber};
% for i = 1:length(data_filtered)
%     data_tmp = data_filtered{i};
%     t_suppress{end+1} = data_tmp.DayNumber;
%     data_tmp = data_tmp{:, 2:end};
%     data_tmp = nthroot(data_tmp, 4);
%     data_suppress{end+1} = data_tmp;
% end

% %% plot suppressed data
% for i = 1:length(data_suppress)
%     figure(i)
%     t_tmp = t_suppress{i};
%     data_tmp = data_suppress{i};
%     for j = 1:num_component
%         plot(t_tmp, data_tmp)
%         hold on
%     end
%     xlim([t_tmp(1), t_tmp(end)])
% end

%% equidistant interpolation
data_total = {};
t_total = {};
for i = 1:length(data_suppress)
    data_tmp = data_suppress{i};
    t_tmp = t_suppress{i};
    data_tmp(isnan(data_tmp)) = 0;
    t_spline = t_tmp(1):1:t_tmp(end);

    % % fill zeros with minimum value(detection limit)
    % cyclopoids_zero_count = sum(data_tmp(:, 1) == 0);
    % cyclopoids_noise = rand(1, cyclopoids_zero_count) * 0.000043;
    % data_tmp(data_tmp(:, 1) == 0, 1) = cyclopoids_noise.';
    % 
    % % fill zeros with minimum value(detection limit)
    % calanoids_zero_count = sum(data_tmp(:, 2) == 0);
    % calanoids_noise = rand(1, calanoids_zero_count) * 0.000043;
    % data_tmp(data_tmp(:, 2) == 0, 2) = calanoids_noise.';
    % 
    % % fill zeros with minimum value(detection limit)
    % rotifers_zero_count = sum(data_tmp(:, 3) == 0);
    % rotifers_noise = rand(1, rotifers_zero_count) * 0.000093;
    % data_tmp(data_tmp(:, 3) == 0, 3) = rotifers_noise.';
    % 
    % % fill zeros with minimum value(detection limit)
    % protozoa_zero_count = sum(data_tmp(:, 4) == 0);
    % protozoa_noise = rand(1, protozoa_zero_count) * 0.000001;
    % data_tmp(data_tmp(:, 4) == 0, 4) = protozoa_noise.';
    % 
    % % fill zeros with minimum value(detection limit)
    % nano_zero_count = sum(data_tmp(:, 5) == 0);
    % nano_noise = rand(1, nano_zero_count) * 0.000243187;
    % data_tmp(data_tmp(:, 5) == 0, 5) = nano_noise.';
    % 
    % % fill zeros with minimum value(detection limit)
    % pico_zero_count = sum(data_tmp(:, 6) == 0);
    % pico_noise = rand(1, pico_zero_count) * 7.56959E-06;
    % data_tmp(data_tmp(:, 6) == 0, 6) = pico_noise.';

    data_total_tmp = zeros(length(t_spline), num_component);
    for j = 1:num_component
        data_spline = interp1(t_tmp, data_tmp(:,j), t_spline.');
        data_total_tmp(:,j) = data_spline;
    end
    data_total{end+1} = data_total_tmp;
    t_total{end+1} = t_spline.';
end

%% plot equidistant data
for i = 1:length(data_total)
    data_before = data_suppress{i};
    t_before = t_suppress{i};
    data_after = data_total{i};
    t_after = t_total{i};
    for j = 1:num_component
        figure(i* 10 + j);
        plot(t_before, data_before(:,j), 'k',t_after, data_after(:,j), 'r')
        legend(labels(j));
    end
end

%% cut time series data
window_size_ori = 120;
overlapping_ratio = 0.1;
% time_interval = 0.5;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

% window_size = (window_size_ori-1) / time_interval + 1;

y_ori = {};
for i = 1:length(data_total)
    start = 1;
    data_tmp = data_total{i};
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

num_fourier = 8;
option = ['fourier', num2str(num_fourier)];

y_fit = {};
t_fit = linspace(0, 1, window_size_ori).';
% t_fit = linspace(0, 1, window_size);

for i = 1:num_data
    y_fit_tmp = y_ori{i};
    for j = 1:num_component
        % y_interp = interp1(t_ori, y_fit_tmp(:,j),t_fit);
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

y_ori_tmp = y_ori{idx};
y_fit_tmp = y_fit{idx};
for i = 1:num_component
    figure(i)
    plot(t_fit, y_ori_tmp(:,i), 'k', t_fit, y_fit_tmp(:,i), 'r-o');
    legend(labels(i))
end

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
y = table2array(data_ori);
y_total = y_fit;
time_interval = 1/window_size_ori;

save('data_cut', 't', 'y', 'y_total', 'time_interval', 'num_data', 'num_component', 'labels')
