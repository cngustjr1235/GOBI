clc;
clear;
close all;

%% load transformed data
load('transformed_data.mat')

%% parameters
t = transformed_data.DayNumber;
data_cut = removevars(transformed_data, 'DayNumber');

%% plot
figure(1)
plot(t, nthroot(data_cut.Rotifers, 4), Color="#D82E32", Marker='o')
xlim([1000, 2000])
legend('Rotifers')

figure(2)
plot(t, nthroot(data_cut.Calanoids, 4), Color="#80307C", Marker='o')
xlim([1000, 2000])
legend('Calanoids')

figure(3)
plot(t, nthroot(data_cut.Picophytoplankton, 4), Color="#2F357E", Marker='o')
xlim([1000, 2000])
legend('Picophytoplankton')

figure(4)
plot(t, nthroot(data_cut.Nanophytoplankton, 4), Color="#009648", Marker='o')
xlim([1000, 2000])
legend('Nanophytoplankton')