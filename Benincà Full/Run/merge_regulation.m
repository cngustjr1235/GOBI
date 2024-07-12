function [regulation_network_1D, regulation_network_2D] = merge_regulation(pair_list_1D, detection_list_1D, pair_list_2D, boot_candidate_list, boot_type_list, surrogate_list, num_component)
%% 1D regulation
regulation_network_1D = zeros(num_component);

for component = 1:num_component
    % indexing all the causes of each target
    index_list = [];
    for pair = 1:length(pair_list_1D(:,1))
        if pair_list_1D(pair,2) == component
            index_list = [index_list, pair];
        end
    end

    % check inferred 1D regulation
    regulation_idx_list = [];
    for pair = index_list
        for type = 1:length(detection_list_1D(1,:))
            if detection_list_1D(pair, type) == 1
                regulation_idx_list = [regulation_idx_list;[pair, type]];
            end
        end
    end

    if ~isempty(regulation_idx_list)
        if length(regulation_idx_list(:,1)) == 1
            pair_idx = regulation_idx_list(1,1);
            type_idx = regulation_idx_list(1,2);

            cause = pair_list_1D(pair_idx,1);
            target = pair_list_1D(pair_idx,2);
            if type_idx == 1
                regulation_network_1D(cause, target) = 1;
            else
                regulation_network_1D(cause, target) = -1;
            end
        end
    end
end

%% 2D regulation  
if isempty(boot_candidate_list)
    num_candidate_boot = 0;
else
    num_candidate_boot = length(boot_candidate_list(:,1)); % candidate for the surrogate test (passing the delta test)
end

% find the index for candidates of 2D regulations
candidate = [];
for pair_idx = 1:num_candidate_boot
    idx1 = find(pair_list_2D(:,1) == boot_candidate_list(pair_idx,1));
    idx2 = find(pair_list_2D(:,2) == boot_candidate_list(pair_idx,2));
    idx3 = find(pair_list_2D(:,3) == boot_candidate_list(pair_idx,3));
    idx = intersect(intersect(idx1, idx2), idx3);
    
    candidate = [candidate; [idx, boot_type_list(pair_idx)]];
end

%% make basis network to find the potential indirect effect (merging every inferred 2D candidates with inferred 1D regulations)
regulation_2D = regulation_network_1D;

for pair_idx = 1:num_candidate_boot
    tar_pair = candidate(pair_idx,1);
    tar_type = candidate(pair_idx,2);
    if tar_type == 1
        type_idx = [1,1];
    elseif tar_type == 2
        type_idx = [1,-1];
    elseif tar_type == 3
        type_idx = [-1,1];
    else
        type_idx = [-1,-1];
    end

    cause1 = pair_list_2D(tar_pair,1);
    cause2 = pair_list_2D(tar_pair,2);
    target  = pair_list_2D(tar_pair,3);

    regulation_2D(cause1, target) = type_idx(1);
    regulation_2D(cause2, target) = type_idx(2);
end

%% check each candidate whether indirect or not
indirect_list = [];
for pair_idx = 1:num_candidate_boot
    indirect_idx = [0,0];
    tar_pair = candidate(pair_idx,1);
    tar_type = candidate(pair_idx,2);
    
    if tar_type == 1
        type_idx = [1,1];
    elseif tar_type == 2
        type_idx = [1,-1];
    elseif tar_type == 3
        type_idx = [-1,1];
    else
        type_idx = [-1,-1];
    end
    cause1 = pair_list_2D(tar_pair,1);
    cause2 = pair_list_2D(tar_pair,2);
    target  = pair_list_2D(tar_pair,3);
    
    % check first causal variable
    result_indirect_1 = [];
    for component = 1:num_component
        if regulation_2D(cause1,component) ~= 0 && component ~= target
            % [ori, prev, curr, tar, type, results, network]
            result_indirect_1 = [result_indirect_1; isIndirect(cause1,cause1,component,target,regulation_2D(cause1, component),result_indirect_1,regulation_2D)];
        end
    end
    if ~isempty(find(result_indirect_1 == type_idx(1)))
        indirect_idx(1) = 1;
    end

    % check second causal variable
    result_indirect_2 = [];
    for component = 1:num_component
        if regulation_2D(cause2, component) ~= 0 && component ~= target
            % [ori, prev, curr, tar, type, results, network]
            result_indirect_2 = [result_indirect_2; isIndirect(cause2,cause2,component,target,regulation_2D(cause2, component),result_indirect_2,regulation_2D)];
        end
    end
    if ~isempty(find(result_indirect_2 == type_idx(2)))
        indirect_idx(2) = 1;
    end
    indirect_list = [indirect_list; indirect_idx];
end

% it is possible to consider all the regulations as potential indirect regulations
%indirect_list = ones(num_candidate_boot, 2);

%% infer the regulation using p_surrogate
% using the results of surrogate test and potential indirect regulations
boot_2D = ones(num_candidate_boot,2);
for pair_idx = 1:num_candidate_boot
    for component = 1:2
        if indirect_list(pair_idx, component) == 1 && surrogate_list(pair_idx, component) > surrogate_list(pair_idx, component+2)
            boot_2D(pair_idx, component) = 0;
        end
    end
end

% merge the results (nan means indirect regulation)
regulation_network_2D = zeros(num_component);
for pair_idx = 1:num_candidate_boot
    cause1 = boot_candidate_list(pair_idx, 1);
    cause2 = boot_candidate_list(pair_idx, 2);
    target = boot_candidate_list(pair_idx, 3);
    if boot_2D(pair_idx, 1) == 1
        if ~isnan(regulation_network_2D(cause1,target))
            if boot_type_list(pair_idx) == 1 || boot_type_list(pair_idx) == 2
                regulation_network_2D(cause1,target) = 1;
            else
                regulation_network_2D(cause1,target) = -1;
            end
        end
    else
        regulation_network_2D(cause1, target) = nan;
    end
    if boot_2D(pair_idx, 2) == 1
        if ~isnan(regulation_network_2D(cause2,target))
            if boot_type_list(pair_idx) == 1 || boot_type_list(pair_idx) == 3
                regulation_network_2D(cause2,target) = 1;
            else
                regulation_network_2D(cause2,target) = -1;
            end
        end
    else
        regulation_network_2D(cause2,target) = nan;
    end
end
for i = 1:num_component
    for j = 1:num_component
        if isnan(regulation_network_2D(i,j))
            regulation_network_2D(i,j) = 0;
        end
    end
end

end

