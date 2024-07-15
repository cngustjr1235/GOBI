clc;
clear;
close all;
addpath('GOBI') 

%% load data
load('data_with_options.mat')

%% parameters
thres_noise = 0; % threshold for regulation-detection function, 0 is default
dimension = 3;   % dimension of the framework

%% all the pairs of 3D regulation (cause1, cause2, cause3, target)

component = [1:num_component];
component_list_dim3_tmp = nchoosek(component, 3);
component_list_dim3 = [];
for i = 1:length(component_list_dim3_tmp(:,1))
    for j = fixed_target
        if ismember(j, component_list_dim3_tmp(i,:))
            if isnan(type_self)
                component_list_dim3 = [component_list_dim3 ; [component_list_dim3_tmp(i,:), j]];
            else
                continue
            end
        else
            component_list_dim3 = [component_list_dim3 ; [component_list_dim3_tmp(i,:), j]];
        end
    end
end

num_pair = length(component_list_dim3(:,1));
num_type = 2.^dimension;

%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
wb = waitbar(0,'Processing');
N = num_data;
Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));

%% from all data, calculate regulation detection score & region
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection score for all data
R_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection region for all data
parfor i = 1:num_data
    send(dq,i)
    y_target = y_total{i};
    t_target = t_total{i};
    
    S_total = zeros(num_pair,num_type); % save regulation-detection score for each data
    R_total = zeros(num_pair,num_type); % save regulation-detection region for each data
    
    for j = 1:length(component_list_dim3(:,1))
        st1 = component_list_dim3(j,1); % index for cause1
        st2 = component_list_dim3(j,2); % index for cause2
        st3 = component_list_dim3(j,3); % index for cause3
        ed = component_list_dim3(j,4);  % index for target
        
        % compute regulation-detection function
        [score_list, t_1, t_2] = RDS_dim3(y_target(:,st1), y_target(:,st2), y_target(:,st3), y_target(:,ed), t_target);
        
        % compute regulation-detection score
        for k = 1:num_type
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            loca_plus = find(score > thres_noise);
            loca_minus = find(score < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s = 1;    
            else
                s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
            end
            r = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            S_total(j,k) = s;
            R_total(j,k) = r;
        end
    end
    S_total_list(:,:,i) = S_total;
    R_total_list(:,:,i) = R_total;    
end

filename = 'RDS_dim3';
save(filename, 'S_total_list', 'R_total_list','component_list_dim3', 'num_component', 'num_pair','num_type','num_data')

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
