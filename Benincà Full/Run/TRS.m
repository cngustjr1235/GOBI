function [TRS_list, detection_list] = TRS(S_total_list, R_total_list, pair_list, thres_S, thres_R, thres_TRS, dimension)
%% compute TRS for each pair of regulation
num_pair = length(pair_list(:,1));
num_type = 2^dimension;
num_data = length(S_total_list(1,1,:));

TRS_list = zeros(num_pair, num_type);
for pair = 1:num_pair
    for type = 1:num_type
        S_tmp = reshape(S_total_list(pair,type,:),[num_data,1]); % import regulation-detection score
        R_tmp = reshape(R_total_list(pair,type,:),[num_data,1]); % import regulation-detection region
        
        S_filtered = S_threshold(S_tmp, thres_S); % test whether S > S^thres
        R_filtered = R_threshold(R_tmp, thres_R); % test whether R < R^thres
        
        if sum(R_filtered) == 0
            TRS_tmp = nan;
        else
            TRS_tmp = sum(S_filtered .* R_filtered) / sum(R_filtered); % compute TRS
        end
        TRS_list(pair, type) = TRS_tmp;
    end
end

%% infer regulation using the criteria TRS > TRS^thres
detection_list = zeros(num_pair, num_type);
for pair = 1:num_pair
    for type = 1:num_type
        if TRS_list(pair,type) >= thres_TRS
            detection_list(pair,type) = 1;
        end
    end
end

end

