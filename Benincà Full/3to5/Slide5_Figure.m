clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
% y = [data_ori.Cyclopoids, data_ori.Rotifers, data_ori.Picophytoplankton];
% labels = ["Cyclopoids", "Rotifers", "Picophytoplankton"];
y = [data_ori.Cyclopoids, data_ori.Picophytoplankton];
labels = ["Cyclopoids", "Picophytoplankton"];
y(isnan(y)) = 0;
t = data_ori.DayNumber;
y = nthroot(y, 4);
num_component = length(y(1,:));

%% equidistant interpolation

t_fit = (t(1):1:t(end)).';
y_fit = zeros(length(t_fit), num_component);
for i = 1:num_component
    y_fit(:,i) = interp1(t, y(:, i), t_fit, 'linear');
end

t_cut1 = t_fit(1:120);
y_cut1 = y_fit(1:120,:);

t_cut2 = t_fit(2539:2658);
y_cut2 = y_fit(2539:2658,:);

for i = 1:num_component
    y_cut1(:,i) = y_cut1(:,i) - min(y_cut1(:,i));
    y_cut1(:,i) = y_cut1(:,i) / max(y_cut1(:,i));

    y_cut2(:,i) = y_cut2(:,i) - min(y_cut2(:,i));
    y_cut2(:,i) = y_cut2(:,i) / max(y_cut2(:,i));
end
y_fit1 = y_cut1;
y_fit2 = y_cut2;
t_fit = t_cut1 / t_cut1(end);

num_fourier = 8;
option = ['fourier', num2str(num_fourier)];

for i = 1:num_component
    fitting = fit(t_fit, y_fit1(:,i), option);
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for k = 1:num_fourier
        fouriers = [fouriers; cos(k*w*t_fit.')];
        fouriers = [fouriers; sin(k*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_fit1(:,i) = coeffs(1:end-1) * fouriers;
    
    fitting = fit(t_fit, y_fit2(:,i), option);
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for k = 1:num_fourier
        fouriers = [fouriers; cos(k*w*t_fit.')];
        fouriers = [fouriers; sin(k*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_fit2(:,i) = coeffs(1:end-1) * fouriers;
end

figure(Position=[0, 0, 500, 400])
plot(t_fit, y_fit1(:,1), LineWidth=2, Color="red", LineStyle="-")
hold on
plot(t_fit, y_fit1(:,2), LineWidth=2, Color="blue", LineStyle="-")
hold on
plot(t_fit, y_fit1(:,3), LineWidth=2, Color="#2A8947", LineStyle="-")

xlim([t_fit(1), t_fit(end)])
ylim([0,1])
fontsize(20, "points")
yticks([0, 1])
xticks([0, 1])
legend(labels)

figure(Position=[0, 0, 500, 400])
plot(t_fit, y_fit2(:,1), LineWidth=2, Color="red", LineStyle="-")
hold on
plot(t_fit, y_fit2(:,2), LineWidth=2, Color="blue", LineStyle="-")
hold on
plot(t_fit, y_fit2(:,3), LineWidth=2, Color="#2A8947", LineStyle="-")

xlim([t_fit(1), t_fit(end)])
ylim([0,1])
fontsize(20, "points")
yticks([0, 1])
xticks([0, 1])
legend(labels)

