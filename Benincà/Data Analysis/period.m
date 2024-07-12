%% To calculate period, use autocorrelation
% (refering to SupplFig8, genetic oscillator

clc;
clear;
close all;

%% import data
load('data_ori.mat')
y = [data_ori.Cyclopoids, data_ori.Rotifers, data_ori.Picophytoplankton];
y(isnan(y)) = 0;
t_tmp = data_ori.DayNumber;

num_component = 3;

%% equidistant interpolation
t_spline = (t_tmp(1):1:t_tmp(end)).';
data_total = zeros(length(t_spline), num_component);
for j = 1:num_component
    data_spline = interp1(t_tmp, y(:,j), t_spline);
    data_total(:,j) = data_spline;
end

y = data_total;
%% for each pair of data, calculate period

length_timeseries = length(y(:,1))-1;
period_list = [];
y_norm = y;
for j = 1:num_component
    y_norm(:,j) = y_norm(:,j) - min(y_norm(:,j));
    y_norm(:,j) = y_norm(:,j) / max(y_norm(:,j));
end

for j = 1:num_component
    [autocor, lag] = xcorr(y_norm(:,j), length_timeseries-1, 'coeff');
    [peak, locs] = findpeaks(autocor);
    if ~isempty(find(locs > length_timeseries))
        tmp = locs(find(locs > length_timeseries));
    else
        tmp = nan;
    end

    per = tmp(1) - length_timeseries;
    period_list = [period_list; per];
end