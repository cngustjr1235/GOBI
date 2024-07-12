function [S_total_list, R_total_list, pair_list] = RDS_2D_new(y_total, t_total, target, num_component, type_self, thres_noise)
%% all the pairs of 2D regulation (cause1, cause2, target)
dimension = 2;   % dimension of the framework
num_data = length(y_total);

component_list = [1:num_component];
cause_pair_list = nchoosek(component_list, 2);
pair_list = [];
for cause_pair = 1:length(cause_pair_list(:,1))
    causes = cause_pair_list(cause_pair, :);
    if ismember(target, causes)
        if isnan(type_self)
            pair_list = [pair_list ; [causes, target]];
        else
            continue
        end
    else
        pair_list = [pair_list ; [causes, target]];
    end
end

num_pair = length(pair_list(:,1));
num_type = 2^dimension;

%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
% wb = waitbar(0,'Processing');
% N = num_data;
% Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));

%% from all data, calculate regulation detection score & region
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection score for all data
R_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection region for all data

parfor data = 1:num_data
    send(dq,data)
    y = cell2mat(y_total(data));
    t = cell2mat(t_total(data));
    
    S_total = zeros(num_pair,num_type); % save regulation-detection score for each data
    R_total = zeros(num_pair,num_type); % save regulation-detection region for each data
    
    for pair = 1:length(pair_list(:,1))
        cause1_index = pair_list(pair,1); % index for cause1
        cause2_index = pair_list(pair,2); % index for cause2
        target_index = pair_list(pair,3); % index for target
        
        % compute regulation-detection function
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim2(y(:,cause1_index), y(:,cause2_index), y(:,target_index), t, nan);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim2(y(:,cause1_index), y(:,cause2_index), y(:,target_index), t, nan);
        else
            [score_list, t_1, t_2] = RDS_dim2(y(:,cause1_index), y(:,cause2_index), y(:,target_index), t, nan);
        end
        
        % compute regulation-detection score
        for type = 1:num_type
            score = reshape(score_list(:,:,type),[length(t),length(t)]);
            loca_plus = find(score > thres_noise);
            loca_minus = find(score < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s = 1;    
            else
                s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
            end
            r = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            S_total(pair,type) = s;
            R_total(pair,type) = r;
        end
    end
    S_total_list(:,:,data) = S_total;
    R_total_list(:,:,data) = R_total;    
end

%% plot heatmap

% plot_RDS_heatmap(dimension, component_list_dim1, S_total_list, labels);

delete(gcp('nocreate'))

%% function for parallel pool
function j = completedJobs(varargin)
    persistent n
    if isempty(n)
        n = 0;
    end
    if numel(varargin) ~=0
    else
        n = n+1;
    end
    j=n;
end


end

