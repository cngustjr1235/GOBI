clc;
clear;
close all;
addpath('../GOBI')

%% load data
load('data_with_options')
load('RDS_dim1')

%% parameters
if thres_L == 0
    thres_L = 0.05;
end
thres_noise = 0;

%% surrogate test
num_boot = 100;
surrogate_list = [];

% disp(i)
st = 1; % causal variable
ed = 3; % target variable
% type_tmp = boot_type_list(i); 
type_tmp = 1;
p_total = [];
% num_data = 50;

%% export data to excel for python surrogate
export_data = zeros(length(t)+1, num_data);
export_data(1,:) = 1:num_data;
for i = 1:num_data
    data_tmp = y_total{i};
    data_tmp(isnan(data_tmp)) = 0;
    cause = data_tmp(:, st);
    export_data(2:end,i) = cause;
end

writematrix(export_data, 'cause_original.xlsx');

%% after surrogate
surrogated = zeros(length(t), num_boot, num_data);
for i = 1:num_data
    disp(i)
    s = readtable('result.xlsx', 'Sheet', ['bin', num2str(i-1)]);
    s_trim = s{2:end, 2:end};
    surrogated(:,:,i) = s_trim;
end

%% surrogate test

s_ori_list = [];
for j = 1:num_data
    y_tmp = cell2mat(y_total(j));    
    C1 = y_tmp(:,st); % import time series of the first variable
    T = y_tmp(:,ed);   % import time series of the target variable
    t_target = t(1:length(y_tmp(:,1)));
    
    boot_tmp = [];
    for k = 1:num_boot
        % bootstrapping of regulation-detection score with shuffled time series of cause 1
        % C1_shuffled = C1(randperm(length(C1))); % shuffle the time seires of the first causal variable
        C1_shuffled = surrogated(:,k,j);

        % compute the regulation-detection score using shuffled time series
        [score_list, t_1, t_2] = RDS_dim1(C1_shuffled, T, t_target, time_interval);
        
        score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
        
        loca_plus = find(score_tmp > thres_noise);
        loca_minus = find(score_tmp < -thres_noise);
        if isempty(loca_plus) && isempty(loca_minus)
            s_tmp_1 = 1;
        else
            s_tmp_1 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
        end
        
        boot_tmp = [boot_tmp ;[s_tmp_1]]; % save the RDS from shuffled time series
    end
    
    % compute the regulation-detection score with original time series
    [score_list, t_1, t_2] = RDS_dim1(C1, T, t_target, time_interval);
    
    score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);

    loca_plus = find(score_tmp > thres_noise);
    loca_minus = find(score_tmp < -thres_noise);
    if isempty(loca_plus) && isempty(loca_minus)
        s_ori = 1;
    else
        s_ori = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
    end
    s_ori_list = [s_ori_list, [s_ori]];
    l_ori = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
 
    % using one sided Z test to compute p value (z score)
    [h1,p1] = ztest(s_ori, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
    p_total = [p_total ; [p1]];
end
if ~isempty(p_total)
    p_tmp_1 = nonzeros(rmmissing(p_total(:,1)));
    
    % combine p-value using Fisher's method
    fisher_tmp_1 = 2* sum(-log(p_tmp_1));
    num_p_1 = length(p_tmp_1);
    
    % combine p_surrogate to define the threshold for the surrogate test
    fisher_thres_1 = chi2cdf(-2*log(p_surrogate)*num_p_1, 2*num_p_1, 'upper'); % threshold of combined p-value for cause1
    fisher_tmp = [chi2cdf(fisher_tmp_1, 2*num_p_1, 'upper'),fisher_thres_1];
else
    fisher_tmp = [0,0,0,0];
end
surrogate_list = [surrogate_list; [fisher_tmp]];

% save the results of surrogate test
filename = ['Surrogate_dim1'];
save(filename,'component_list_dim1','num_data','surrogate_list')

