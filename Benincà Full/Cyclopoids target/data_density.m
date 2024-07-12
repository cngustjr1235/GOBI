clc;
clear;
close all;
load('data_cut.mat')

%% load raw data
load('data_ori.mat')
t = data_ori.DayNumber;
t = t(1:end, :);

%% data density

% must be same with Step0_2_Original
window_size_ori = 120;
overlapping_ratio = 0.9;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

density = [];
start = 1;
length_timeseries = t(end);
while(1)
    if start + window_size_ori > length_timeseries
        start = length_timeseries - window_size_ori + 1;
    end
    data_cnt = sum(t >= start & t < start + window_size_ori);
    density = [density, [data_cnt / window_size_ori]];
    if start + window_size_ori > length_timeseries
        break;
    end
    start = start + window_move_ori;
end


