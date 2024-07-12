clc;
clear;
close all;

%% load raw data
load('original_data.mat')

%% remove duplicate line?
data(701,:) = [];

%% parameters
t = data.DayNumber;
data_cut = removevars(data, 'DayNumber');
num_species = width(data_cut) - 1;

%% plot 11 species
figure(1)
for i = 1:num_species
    plot(t, nthroot(data_cut{:, i}, 4))
    hold on
end
xlim([0, t(end)])
ylim([0, 6.5])
xlabel('Time (days)')
ylabel('Species')
legend(data_cut.Properties.VariableNames)

%% plot 4 species
figure(2)
plot(t, nthroot(data_cut.Rotifers, 4), Color="#D82E32", LineWidth=1)
hold on
plot(t, nthroot(data_cut.CalanoidCopepods, 4), Color='#80307C', LineWidth=1)
hold on
plot(t, nthroot(data_cut.Picophytoplankton, 4), Color='#2F357E', LineWidth=1)
hold on
plot(t, nthroot(data_cut.Nanophytoplankton, 4), Color='#009648', LineWidth=1)

xlim([1550, 1750])
ylim([0, 3.0])
xlabel('Time (days)')
ylabel('Plankton')
legend('Rotifers', 'CalanoidCopepods', 'Picophytoplankton', 'Nanophytoplankton')

%% plot splined data
load('transformed_data.mat')

figure(3)
plot(t, nthroot(data_cut.Rotifers, 4),'ob')
hold on

splined_t = transformed_data.DayNumber;
plot(splined_t, transformed_data.Rotifers, Color='#D82E32', LineStyle='-', Marker='.', LineWidth=1)
xlim([1550,1750])
xlabel('Time (days)')
ylabel('Rotifers')
legend('Raw data', 'Interpolated')