clc;
clear;
close all;

%% load raw data
load('data_ori.mat')
t = data_ori.DayNumber;
y_raw = horzcat(data_ori.Cyclopoids, data_ori.Rotifers, data_ori.Protozoa, data_ori.Picophytoplankton);
y_raw(isnan(y_raw)) = 0;
y = nthroot(y_raw, 4);

%% parameters 1. interpolation

% interpolation method
% method = 1: linear interpolation
% method = 2: spline interpolation
% method = 3: fourier interpolation, 
%             In this case, users have to choose the order of fourier method (1~8)

method = 3;
num_fourier = 8;

%% parameters 2. cut the time series data
num_component = length(y(1,:));        % number of component in this system
window_size_ori = 140;     % For oscillatory data, 1 period is recommended
overlapping_ratio = 0.9;  % overlapping ratio of moving window technique

% choose sampling rate for interpolation. 
% window_size_ori/time_interval becomes number of time points per cutted data
% window_size_ori/time_interval = 100 is recommended
% window_size_ori/time_interval is high (low) make the inference accurate (less accurate) and slow (fast)
time_interval = 1/2;       

%% process parameters
window_size = (window_size_ori - 1) / time_interval + 1;
window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));
window_move = (window_move_ori - 1) / time_interval + 1;


%% plot raw data
figure(1)
for i = 1:num_component
    plot(t, y(:,i))
    hold on
end
xlim([0,t(end)])
ylim([0,max(max(y))])
legend('Rotifers', 'Calanoids', 'Picophytoplankton', 'Nanophytoplankton')

%% cut data
disp('cut data...')
y_total_ori = {};
start = 1;
length_timeseries = length(y(:,1));
while(1)
    if start + window_size_ori > length_timeseries
        break
    end
    y_tmp = y(start:start+window_size_ori-1, :);
    for i = 1:num_component
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i));
    end
    y_total_ori{end+1} = y_tmp;
    start = start + window_move_ori;
end
num_data = length(y_total_ori);

%% plot cutted data
t = linspace(0, 1, window_size_ori).';
% for i = 1:num_data
%     figure(i + 1);
%     y_tmp = y_total_ori{i};
%     for j = 1:num_component
%         plot(t, y_tmp(:,j));
%         hold on;
%     end
%     xlim([0,t(end)])
%     ylim([0,max(max(y_tmp))])
%     legend('Rotifers', 'Calanoids', 'Picophytoplankton', 'Nanophytoplankton')
% end

%% interpolate data

y_total = {};

t_fit = linspace(0, 1, window_size).';

for i = 1:num_data
    y_fit = zeros(length(t_fit), num_component);
    y_tmp = y_total_ori{i};
    if method == 1 % linear
        for j = 1:num_component
            y_fit(:,j) = interp1(t, y_tmp(:,j),t_fit,'linear');
        end
    elseif method == 2 % spline
        for j = 1:num_component
            y_fit(:,j) = interp1(t, y_tmp(:,j),t_fit,'spline');
        end
    elseif method == 3 % fourier
        for j = 1:num_component
            y_fit(:,j) = interp1(t, y_tmp(:,j),t_fit,'linear');
        end
        option = ['fourier',num2str(num_fourier)];
        for j = 1:num_component
            fitting = fit(t_fit,y_fit(:,j),option);
            w = fitting.w;
            fouriers = ones(1,length(t_fit));
            for k = 1:num_fourier
                fouriers = [fouriers; cos(k*w*t_fit.')];
                fouriers = [fouriers; sin(k*w*t_fit.')];
            end
            coeffs = coeffvalues(fitting);
            y_fit(:,j) = coeffs(1:end-1) * fouriers;
        end
    end
    y_fit(y_fit < 0) = 0;
    y_total{end+1} = y_fit;
end

%% plot interpolated data
% for i = 1:num_data
%     y_tmp_ori = y_total_ori{i};
%     y_tmp = y_total{i};
%     for j = 1:num_component
%         figure(i * 10 + j);
%         plot(t, y_tmp_ori(:,j), t_fit, y_tmp(:,j))
%         xlim([0,t_fit(end)])
%         ylim([0,1])
%         if j == 1
%             title('Rotifer ',i)
%         elseif j == 2
%             title('Calanoids ', i);
%         elseif j == 3
%             title('Pico ', i);
%         else 
%             title('Nano ', i);
%         end
%     end
%     % legend('Rotifers', 'Calanoids', 'Picophytoplankton', 'Nanophytoplankton')
% end

%% Save data
time_interval = 1/window_size;
save('data_cut','t','y','y_total','time_interval','num_data','num_component')