function [delta_candidate_list, delta_type_list, delta_list] = delta_test_2D(S_total_list, R_total_list, detection_list, pair_list, thres_R)
%% find candidate for delta test
delta_candidate_list = [];
delta_type_list = [];

dimension = 2;
num_pair = length(pair_list(:,1));
num_type = 2^dimension;
num_data = length(S_total_list(1,1,:));

for pair = 1:num_pair
    for type = 1:num_type
        if detection_list(pair, type) == 1
            delta_candidate_list = [delta_candidate_list ; [pair_list(pair,:)]];
            delta_type_list = [delta_type_list ; type];
        end
    end
end

if isempty(delta_candidate_list)
    num_candidate_delta = 0;
else
    num_candidate_delta = length(delta_candidate_list(:,1));
end

%% delta test
delta_list = zeros(num_candidate_delta, 2); % update the result of delta

% correspoding type used for delta test
% ex1> to test type (+,+), delta test use scores of (-,+) and (+,-)
%      thus, the first row of cor_type is [3,2]
% ex2> to test type (+,-), delta test use scores of (-,-) and (+,+)
%      thus, the first row of cor_type is [4,1]
% ex3> to test type (-,+), delta test use scores of (+,+) and (-,-)
%      thus, the first row of cor_type is [1,4]
% ex4> to test type (-,-), delta test use scores of (+,-) and (-,+)
%      thus, the first row of cor_type is [2,3]

cor_type = [
    [3,2];
    [4,1];
    [1,4];
    [2,3]];

if num_candidate_delta == 0
    disp('no candidate for delta')
else
for candidate_pair = 1:num_candidate_delta
    % find index for pair
    cause1 = delta_candidate_list(candidate_pair,1);
    cause2 = delta_candidate_list(candidate_pair,2);
    target = delta_candidate_list(candidate_pair,3);
    idx_cause1 = find(pair_list(:,1) == cause1);
    idx_cause2 = find(pair_list(:,2) == cause2);
    idx_target = find(pair_list(:,3) == target);
    idx = intersect(intersect(idx_cause1, idx_cause2),idx_target);
    
    % import regulation-detection score and region
    type_tmp = delta_type_list(candidate_pair);
    S_tmp_ori = reshape(S_total_list(idx,type_tmp,:), [num_data,1]);
    S_tmp_1 = reshape(S_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    S_tmp_2 = reshape(S_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    R_ori = reshape(R_total_list(idx,type_tmp,:), [num_data,1]);
    R_tmp_1 = reshape(R_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    R_tmp_2 = reshape(R_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    
    % use only when R > R^thres
    R_processed_ori = R_threshold(R_ori,thres_R);
    R_processed_1 = R_threshold(R_tmp_1,thres_R);
    R_processed_2 = R_threshold(R_tmp_2,thres_R);
    S_processed_ori = S_tmp_ori .* R_processed_ori;
    S_processed_1 = S_tmp_1 .* R_processed_1;
    S_processed_2 = S_tmp_2 .* R_processed_2;
    
    S_processed_1(find(S_processed_1 == 0)) = NaN;
    S_processed_2(find(S_processed_2 == 0)) = NaN;
    S_processed_ori(find(S_processed_ori == 0)) = NaN;
    
    % compute regulation-delta function
    delta_1 = S_processed_ori - S_processed_1;
    delta_2 = S_processed_ori - S_processed_2;
    
    p1 = nan;
    p2 = nan;
    
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
            [p2,h2,stats2] = signrank(rmmissing(delta_2),-1e-3,'tail','right');
        else
            p2 = length(find(rmmissing(delta_2) < -1e-2)) / length(rmmissing(delta_2));
        end
    end
    delta_list(candidate_pair,:) = [p1,p2];
end

end

