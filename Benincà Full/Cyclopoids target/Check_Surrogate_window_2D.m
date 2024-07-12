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
thres_noise = 0;
dimension = 2;

%% Pro, Rot -> Cyc only

S_target = reshape(S_total_list(3,1,:), [1, num_data]);
L_target = reshape(L_total_list(3,1,:), [1, num_data]);

S_target(L_target < thres_L) = nan; % When region is zero

%% surrogate test
num_boot = 1;
surrogate_list = [];

st1 = component_list_dim2(3,1); % first causal variable
st2 = component_list_dim2(3,2); % second causal variable
ed = component_list_dim2(3,3);  % target variable
type_tmp = 1;   
p_total = [];
p_result_total = [];
% for j = 1:num_data
for j = 52:52
    y_tmp = cell2mat(y_total(j));    
    C1 = y_tmp(:,st1); % import time series of the first variable
    C2 = y_tmp(:,st2); % import time series of the second variable
    T = y_tmp(:,ed);   % import time series of the target variable
    t_target = t(1:length(y_tmp(:,1)));
    
    boot_tmp = [];
    for k = 1:num_boot
        % bootstrapping of regulation-detection score with shuffled time series of cause 1
        C1_shuffled = C1(randperm(length(C1))); % shuffle the time seires of the first causal variable
        
        % compute the regulation-detection score using shuffled time series
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim2(C1_shuffled, C2, T, t_target, time_interval);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim2(C1_shuffled, C2, T, t_target, time_interval);
        else
            [score_list, t_1, t_2] = RDS_dim2(C1_shuffled, C2, T, t_target, time_interval);
        end
        
        score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
        
        loca_plus = find(score_tmp > thres_noise);
        loca_minus = find(score_tmp < -thres_noise);
        if isempty(loca_plus) && isempty(loca_minus)
            s_tmp_1 = 1;
        else
            s_tmp_1 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
        end
        
        % bootstrapping of regulation-detection score with shuffled time series of cause 1
        C2_shuffled = C2(randperm(length(C2))); % shuffle the time seires of the second causal variable
        
        % compute the regulation-detection score using shuffled time series
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2_shuffled, T, t_target, time_interval);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim2(C1, C2_shuffled, T, t_target, time_interval);
        else
            [score_list, t_1, t_2] = RDS_dim2(C1, C2_shuffled, T, t_target, time_interval);
        end
        
        score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
        
        loca_plus = find(score_tmp > thres_noise);
        loca_minus = find(score_tmp < -thres_noise);
        if isempty(loca_plus) && isempty(loca_minus)
            s_tmp_2 = 1;
        else
            s_tmp_2 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
        end
        
        boot_tmp = [boot_tmp ;[s_tmp_1,s_tmp_2]]; % save the RDS from shuffled time series
    end
    
    % compute the regulation-detection score with original time series
    if type_self == -1
        [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2, T, t_target, time_interval);
    elseif type_self == 1
        [score_list, t_1, t_2] = RDS_ps_dim2(C1, C2, T, t_target, time_interval);
    else
        [score_list, t_1, t_2] = RDS_dim2(C1, C2, T, t_target, time_interval);
    end
    
    score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);

    loca_plus = find(score_tmp > thres_noise);
    loca_minus = find(score_tmp < -thres_noise);
    if isempty(loca_plus) && isempty(loca_minus)
        s_ori = 1;
    else
        s_ori = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
    end
    l_ori = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
 
    % using one sided Z test to compute p value (z score)
    [h1,p1] = ztest(s_ori, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
    [h2,p2] = ztest(s_ori, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
    p_total = [p_total ; [p1,p2]];
    if p1 < p_surrogate && p2 < p_surrogate
        p_result_total = [p_result_total ; [1]];
    else
        p_result_total = [p_result_total ; [0]];
    end
end

%% Draw heatmap

c_nan = [0.5 0.5 0.5];
font_s = 14;
tmp = 1 - linspace(0, 1, 251)';
cmap_score = vertcat([[ones(251,1)], [1-tmp], [1-tmp]], [[tmp],[tmp],[ones(251,1)]]);


figure(Position=[0,0,1000,80]);
S_target(p_result_total == 0) = nan;
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