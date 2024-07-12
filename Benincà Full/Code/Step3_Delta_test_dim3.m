clc;
clear;
close all;
addpath('GOBI')

%% load data
load('data_with_options')
load('RDS_dim3')

%% parameters
if thres_L == 0
    thres_L = 0.01;
end

%% compute TRS
TRS_total = zeros(num_pair, num_type);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total_list(i,j,:),[num_data,1]);
        L_tmp = reshape(L_total_list(i,j,:),[num_data,1]);
        
        S_processed = S_threshold(S_tmp, thres_S);
        L_processed = L_threshold(L_tmp, thres_L);
        
        if sum(L_processed) == 0
            TRS_tmp = 0;
        else
            TRS_tmp = sum(S_processed .* L_processed) / sum(L_processed);
        end
        TRS_total(i,j) = TRS_tmp;
    end
end

%% find candidate for delta test
delta_candidate_list = [];
delta_type_list = [];
for i = 1:num_pair
    for k = 1:num_type
        if TRS_total (i,k) > thres_TRS
            delta_candidate_list = [delta_candidate_list ; [component_list_dim3(i,:)]];
            delta_type_list = [delta_type_list ; k];
        end
    end
end
num_candidate_delta = length(delta_candidate_list(:,1));

%% delta test
delta_list = zeros(num_candidate_delta, 3); % update the result of delta

% corresponding type used for delta test
% ex1> to test type (+,+,+), delta test use scores of (-,+,+) and (+,-,+)
%      and (+,+,-) thus, the cor_type is [5,3,2]
% ex2> to test type (+,+,-), delta test use scores of (-,+,-) and (+,-,-)
%      and (+,+,+) thus, the cor_type is [6,4,1]
% ex3> to test type (+,-,+), delta test use scores of (-,-,+) and (+,+,+)
%      and (+,-,-) thus, the cor_type is [7,1,4]
% ex4> to test type (+,-,-), delta test use scores of (-,-,-) and (+,+,-)
%      and (+,-,+) thus, the cor_type is [8,2,3]
% ex5> to test type (-,+,+), delta test use scores of (+,+,+) and (-,-,+)
%      and (-,+,-) thus, the cor_type is [1,7,6]
% ex6> to test type (-,+,-), delta test use scores of (+,+,-) and (-,-,-)
%      and (-,+,+) thus, the cor_type is [2,8,5]
% ex7> to test type (-,-,+), delta test use scores of (+,-,+) and (-,+,+)
%      and (-,-,-) thus, the cor_type is [3,5,8]
% ex8> to test type (-,-,-), delta test use scores of (+,-,-) and (-,+,-)
%      and (-,-,+) thus, the cor_type is [4,6,7]

cor_type = [
    [5,3,2];
    [6,4,1];
    [7,1,4];
    [8,2,3];
    [1,7,6];
    [2,8,5];
    [3,5,8];
    [4,6,7]];

component_list = component_list_dim3;

if num_candidate_delta == 0
    disp('no candidate for delta')
else
for i = 1:num_candidate_delta
    % find index for pair
    st1 = delta_candidate_list(i,1);
    st2 = delta_candidate_list(i,2);
    st3 = delta_candidate_list(i,3);
    ed = delta_candidate_list(i,4);
    idx_st1 = find(component_list(:,1) == st1);
    idx_st2 = find(component_list(:,2) == st2);
    idx_st3 = find(component_list(:,3) == st3);
    idx_ed = find(component_list(:,4) == ed);
    idx = intersect(intersect(intersect(idx_st1, idx_st2), idx_st3),idx_ed);
    
    % import regulation-detection score and region
    type_tmp = delta_type_list(i);
    S_tmp_ori = reshape(S_total_list(idx,type_tmp,:), [num_data,1]);
    S_tmp_1 = reshape(S_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    S_tmp_2 = reshape(S_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    S_tmp_3 = reshape(S_total_list(idx,cor_type(type_tmp,3),:), [num_data,1]);
    L_ori = reshape(L_total_list(idx,type_tmp,:), [num_data,1]);
    L_tmp_1 = reshape(L_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    L_tmp_2 = reshape(L_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    L_tmp_3 = reshape(L_total_list(idx,cor_type(type_tmp,3),:), [num_data,1]);
    
    % use only when R > R^thres
    L_processed_ori = L_threshold(L_ori,thres_L);
    L_processed_1 = L_threshold(L_tmp_1,thres_L);
    L_processed_2 = L_threshold(L_tmp_2,thres_L);
    L_processed_3 = L_threshold(L_tmp_3,thres_L);
    S_processed_ori = S_tmp_ori .* L_processed_ori;
    S_processed_1 = S_tmp_1 .* L_processed_1;
    S_processed_2 = S_tmp_2 .* L_processed_2;
    S_processed_3 = S_tmp_3 .* L_processed_3;

    S_processed_1(find(S_processed_1 == 0)) = NaN;
    S_processed_2(find(S_processed_2 == 0)) = NaN;
    S_processed_3(find(S_processed_3 == 0)) = NaN;
    S_processed_ori(find(S_processed_ori == 0)) = NaN;
    
    % compute regulation-delta function
    delta_1 = S_processed_ori - S_processed_1;
    delta_2 = S_processed_ori - S_processed_2;
    delta_3 = S_processed_ori - S_processed_3;
    
    p1 = nan;
    p2 = nan;
    p3 = nan;
    
    % if number of data < 25, we compute the ratio of delta >= 0
    % if number of data > 25, we perform Wilcoxon signed rank test
    
    if ~isempty(rmmissing(delta_1))
        if length(rmmissing(delta_1)) >= 25
            [p1,h1,stats1] = signrank(rmmissing(delta_1),-1e-2,'tail','right');
        else
            p1 = length(find(rmmissing(delta_1) < -1e-2)) / length(rmmissing(delta_1));
        end
    end
    if ~isempty(rmmissing(delta_2))
        if length(rmmissing(delta_2)) >= 25
            [p2,h2,stats2] = signrank(rmmissing(delta_2),-1e-2,'tail','right');
        else
            p2 = length(find(rmmissing(delta_2) < -1e-2)) / length(rmmissing(delta_2));
        end
    end
    if ~isempty(rmmissing(delta_3))
        if length(rmmissing(delta_3)) >= 25
            [p3,h3,stats3] = signrank(rmmissing(delta_3),-1e-2,'tail','right');
        else
            p3 = length(find(rmmissing(delta_3) < -1e-2)) / length(rmmissing(delta_3));
        end
    end
    delta_list(i,:) = [p1,p2,p3];
end
end
save('Delta_dim3', 'delta_list', 'delta_candidate_list', 'delta_type_list', 'TRS_total', 'num_candidate_delta')

