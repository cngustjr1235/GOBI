clc;
clear;
close all;
addpath('../')

%% load raw data
load('data_ori.mat')
y  = [data_ori.Cyclopoids, data_ori.Rotifers, data_ori.Picophytoplankton];
labels = ["Cyclopoids", "Rotifers", "Picophytoplankton"];
num_component = length(labels);

%% suppress data by 1/4
data_suppress = y;
t = data_ori.DayNumber;
data_suppress = nthroot(data_suppress, 4);
data_suppress(isnan(data_suppress)) = 0;


%% equidistant interpolation
t_fit = (t(1):1:t(end)).';
data_fit = zeros(length(t_fit), num_component);
for i = 1:num_component
    data_fit(:,i) = interp1(t, data_suppress(:, i), t_fit, 'linear');
end

%% cut time series data
window_size_ori = 120;
overlapping_ratio = 0.9;

window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));

y_ori = {};
start = 1;
data_tmp = data_fit;
length_timeseries = length(data_tmp(:,1));
while(1)
    if start + window_size_ori > length_timeseries
        start = length_timeseries - window_size_ori + 1;
    end
    y_tmp = data_tmp(start:start+window_size_ori-1, :);
    for j = 1:num_component
        % normalize data
        y_tmp(:,j) = y_tmp(:,j) - min(y_tmp(:,j));
        if max(y_tmp(:,j)) ~= 0
            y_tmp(:,j) = y_tmp(:,j) / max(y_tmp(:,j));
        end
    end
    y_ori{end+1} = y_tmp;
    if start + window_size_ori > length_timeseries
        break;
    end
    start = start + window_move_ori;
end

num_data = length(y_ori);

%% interpolate data

num_fourier = 8;
option = ['fourier', num2str(num_fourier)];

y_fit = {};
t_fit = linspace(0, 1, window_size_ori).';

for i = 1:num_data
    y_fit_tmp = y_ori{i};
    for j = 1:num_component
        fitting = fit(t_fit, y_fit_tmp(:, j), option);
        w = fitting.w;
        fouriers = ones(1, length(t_fit));
        for k = 1:num_fourier
            fouriers = [fouriers; cos(k*w*t_fit.')];
            fouriers = [fouriers; sin(k*w*t_fit.')];
        end
        coeffs = coeffvalues(fitting);
        y_fit_tmp(:,j) = coeffs(1:end-1) * fouriers;
    end
    y_fit{end+1} = y_fit_tmp;
end

%% options

% defaults
noise_level = 30; % compute noise level using residuals
thres_L = 0;     % threshold for regulation-detection region
                 % thres_L = 0 represents default (decreases as dimension increases
% thres_S = 0.9 - 0.005 * noise_level; % threshold for regulation-detection score
thres_S = 0.5;

type_self = nan;
max_D = num_component - 1;
time_interval = 1/window_size_ori;


%% parameters
thres_noise = 0; % threshold for regulation-detection function, 0 is default
dimension = 1;   % dimension of the framework

%% all the pairs of 1D regulation (cause, target)

component = [1:num_component];
component_list_dim1_tmp = nchoosek(component, 1);
component_list_dim1 = [];
for i = 1:length(component_list_dim1_tmp)
    for j = 1:num_component
        if i == j
            if isnan(type_self)
                component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i), j]];
            else
                continue
            end
        else
            component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i), j]];
        end
    end
end

num_pair = length(component_list_dim1(:,1));
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
L_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection region for all data
parfor i = 1:num_data
    send(dq,i)
    y_target = cell2mat(y_fit(i));
    t_target = t;
    
    S_total = zeros(num_pair,num_type); % save regulation-detection score for each data
    L_total = zeros(num_pair,num_type); % save regulation-detection region for eacg data
    
    for j = 1:length(component_list_dim1(:,1))
        st = component_list_dim1(j,1); % index for cause
        ed = component_list_dim1(j,2); % index for target
        
        % compute regulation-detection function
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        else
            [score_list, t_1, t_2] = RDS_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        end
        
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
            l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            S_total(j,k) = s;
            L_total(j,k) = l;
        end
    end
    S_total_list(:,:,i) = S_total;
    L_total_list(:,:,i) = L_total;    
end

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


