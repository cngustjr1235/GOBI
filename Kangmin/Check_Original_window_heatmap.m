clc;
clear;
close all;
addpath('GOBI')
addpath('../Utility')

%% load data
load('data_with_options')
load('RDS_dim2')

%% parameter
if thres_L == 0
    thres_L = 0.05;
end
dimension = 1;

%% Pro, Rot -> Cyc only

S_target = reshape(S_total_list(3,1,:), [1, num_data]);
L_target = reshape(L_total_list(3,1,:), [1, num_data]);

S_target(L_target < thres_L) = nan; % When region is zero

%% Draw heatmap



c_nan = [0.5 0.5 0.5];
font_s = 14;
tmp = 1 - linspace(0, 1, 251)';
cmap_score = vertcat([[ones(251,1)], [1-tmp], [1-tmp]], [[tmp],[tmp],[ones(251,1)]]);


figure(Position=[0,0,1000,80]);
h = heatmap(S_target);

score_range = [-1,1];
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = score_range;
h.FontSize = font_s;
h.MissingDataColor = c_nan;
% h.XData = x_var;
% h.YData = y_var;
% h.YLabel = y_label;
h.XDisplayLabels = repmat({''}, 1, length(h.XDisplayData));
h.YDisplayLabels = ' ';
s = struct(h);
cbh = s.Colorbar;
set(cbh, 'YTick', score_range)