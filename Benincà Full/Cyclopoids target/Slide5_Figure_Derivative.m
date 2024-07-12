clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
% y = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
y = [data_ori.Cyclopoids,data_ori.Rotifers];

labels = ["Derivative of Cyclopoids", "Rotifers"];
y(isnan(y)) = 0;
t = data_ori.DayNumber;
% y = nthroot(y, 4);
num_component = length(y(1,:));

%% equidistant interpolation

t_fit = (t(1):1:t(end)).';
y_fit = zeros(length(t_fit), num_component);
for i = 1:num_component
    y_fit(:,i) = interp1(t, y(:, i), t_fit, 'linear');
end

t_cut1 = t_fit(490:600);
y_cut1 = y_fit(490:600,:);

for i = 1:num_component
    y_cut1(:,i) = y_cut1(:,i) - min(y_cut1(:,i));
    if max(y_cut1(:,i)) ~= 0
        y_cut1(:,i) = y_cut1(:,i) / max(y_cut1(:,i));
    end
end
y_fit1 = y_cut1;

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
end

cycdot = gradient(y_fit1(:,1), t_fit(2)-t_fit(1));

t_fit = t_fit(1:100);
cycdot = cycdot(1:100);

t_fit = t_fit/max(t_fit);

figure(Position=[0, 0, 1000, 400])
% colororder([255 0 0; 0 0 0]/255);
% yyaxis left
plot(t_fit, cycdot, LineWidth=2, Color=[0 122 255]/255, LineStyle="-")
% hold on
% yyaxis right
% plot(t_fit, y_fit1(:,2), LineWidth=2, Color="#2A8947", LineStyle="-")
% plot(t_fit, y_fit1(:,2), LineWidth=2, Color="blue", LineStyle="-")

xlim([t_fit(1), t_fit(end)])
% ylim([0,1])
fontsize(20, "points")

xticks([0, 1])
yticks([-60, -20, 20, 60])
% legend(labels)

%% compute regulation-detection function
t0 = t_fit;
t_target = t_fit;
time_interval = t_fit(2) - t_fit(1);
[score_list, t_1, t_2] = RDS_dim1(y_fit(:,2), y_fit(:,1), t_target, time_interval);
score_target = reshape(score_list(:,:,2), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1c
c_nan = [0.9 0.9 0.9];
font_s = 14;
tmp = linspace(0, 1, 501)';
cmap_score = [[ones(500,1);1-tmp],[tmp(1:end-1);1-tmp],[tmp; ones(500,1)]];
cmap_value = [[ones(500,1);1-tmp],[tmp*0.3+0.7; (1-tmp(1:end-1))*0.6+0.4],[tmp(1:end-1);1-tmp]];

figure(3)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, 1) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1c')

