clc;
clear;
close all;
addpath('GOBI')

%% load data
load('data_with_options')
load('RDS_dim3')

%% parameter
if thres_L == 0
    thres_L = 0.01;
end
dimension = 3;
%% compute TRS for each pair of 3D regulation
TRS_total = zeros(num_pair, num_type);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total_list(i,j,1:end),[num_data,1]); % import regulation-detection score
        L_tmp = reshape(L_total_list(i,j,1:end),[num_data,1]); % import regulation-detection region
        
        S_processed = S_threshold(S_tmp, thres_S); % test whether S > S^thres
        L_processed = L_threshold(L_tmp, thres_L); % test whether R < R_thres
        
        if sum(L_processed) == 0
            TRS_tmp = nan;
        else
            TRS_tmp = sum(S_processed .* L_processed) / sum(L_processed); % compute TRS
        end
        TRS_total(i,j) = TRS_tmp;
    end
end

%% infer 3D regulation using the criteria TRS > TRS^thres
regulation_3dim = zeros(num_pair,2^dimension);
for i = 1:num_pair
    for j = 1:num_type
        if TRS_total(i,j) >= thres_TRS
            regulation_3dim(i,j) = 1;
        end
    end
end

%% save TRS & inferred regulation
TRS_total_dim3 = TRS_total;

filename = ['TRS_dim3'];
save(filename, 'regulation_3dim','TRS_total_dim3','component_list_dim3')

%% plot heatmap

c_nan = [0.9 0.9 0.9];
font_s = 14;
tmp = 1 - linspace(0, 1, 501)';
cmap_score = [[tmp],[tmp],[ones(501,1)]];

x_var = ["+ + +","+ + -", "+ - +", "+ - -", "- + +", "- + -", "- - +", "- - -"];
y_var = [];
for i = 1:length(TRS_total(:,1))
    y_label = [];
    for j = 1:dimension+1
        idx = component_list_dim3(i, j);
        y_label = [y_label, labels(idx)];
    end
    y_var = [y_var, strjoin(y_label, ', ')];
end

h = heatmap(TRS_total);

TRS_range = [0,1];
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = TRS_range;
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XData = x_var;
h.YData = y_var;
s = struct(h);
cbh = s.Colorbar;
set(cbh, 'YTick', TRS_range)