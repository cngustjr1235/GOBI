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
dimension = 2;
%% compute TRS for each pair of 1D regulation
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

%% infer 2D regulation using the criteria TRS > TRS^thres
regulation_2dim = zeros(num_pair,2^dimension);
for i = 1:num_pair
    for j = 1:num_type
        if TRS_total(i,j) >= thres_TRS
            regulation_2dim(i,j) = 1;
        end
    end
end

%% save TRS & inferred regulation
TRS_total_dim2 = TRS_total;

filename = ['TRS_dim2'];
save(filename, 'regulation_2dim','TRS_total_dim2','component_list_dim2')

%% plot heatmap

plot_heatmap(dimension, component_list_dim2, TRS_total, labels);

[Ssorted, I] = sortrows(component_list_dim2, 3);
TRS_sort_total = TRS_total;
for i = 1:2^dimension
    col = TRS_total(:, i);
    TRS_sort_total(:,i) = col(I);
end
plot_heatmap(dimension, Ssorted, TRS_sort_total, labels);