clc;
clear;
close all;

load('power_series_area-pfc_mouse-4_session-1.mat');

cycle = 1;
start = 1;

y_1 = y(:,1);
y_rms = [];
t_rms = [];
length_timeseries = length(y_1);
while(1)
    if start > length_timeseries
        break;
    end
    y_range = y_1(start:end);
    t_range = t(start:end);
    t_bool = t_range < cycle * 0.15;
    t_cnt = sum(t_bool);
    y_range_only = y_range(1:t_cnt);
    y_rms = [y_rms, [rms(y_range_only)]];
    t_rms = [t_rms, [0.15 * (cycle-1)]];
    start = start + t_cnt;
    cycle = cycle + 1;
end

plot(t_rms, y_rms);