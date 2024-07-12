clc;
clear;
close all;
addpath('../GOBI') 

%% load data
load('data_with_options.mat')

%% parameters
thres_noise = 0; % threshold for regulation-detection function, 0 is default
dimension = 2;   % dimension of the framework
range = max(t);

%% parameter for figure
c_nan = [0.9 0.9 0.9];
font_s = 14;
tmp = linspace(0, 1, 501)';
cmap_score = [[ones(500,1);1-tmp],[tmp(1:end-1);1-tmp],[tmp; ones(500,1)]];
% cmap_score = [[ones(500,1) * 0.2], [ones(500,1) * 0.2], [ones(500,1) * 0.2]];
cmap_value = [[ones(500,1);1-tmp],[tmp*0.3+0.7; (1-tmp(1:end-1))*0.6+0.4],[tmp(1:end-1);1-tmp]];

%% all the pairs of 2D regulation (cause1, cause2, target)

component = [1:num_component];
component_list_dim2_tmp = nchoosek(component, 2);
component_list_dim2 = [];
for i = 1:length(component_list_dim2_tmp(:,1))
    for j = 1:1
        if ismember(j, component_list_dim2_tmp(i,:))
            if isnan(type_self)
                component_list_dim2 = [component_list_dim2 ; [component_list_dim2_tmp(i,:), j]];
            else
                continue
            end
        else
            component_list_dim2 = [component_list_dim2 ; [component_list_dim2_tmp(i,:), j]];
        end
    end
end

num_pair = length(component_list_dim2(:,1));
num_type = 2.^dimension;

%% from all data, calculate regulation detection score & region
disp('calculate regulation detection score...')

data_index = 21;
component_index = 3;
type_index = 1; % type index

y_target = cell2mat(y_total(data_index));
t_target = t;

st1 = component_list_dim2(component_index,1); % index for cause1
st2 = component_list_dim2(component_index,2); % index for cause2
ed = component_list_dim2(component_index,3);  % index for target

% compute regulation-detection function
if type_self == -1
    [score_list, t_1, t_2] = RDS_ns_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
elseif type_self == 1
    [score_list, t_1, t_2] = RDS_ps_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
else
    [score_list, t_1, t_2] = RDS_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
end

% plot regulation-detection function
score = reshape(score_list(:,:,type_index),[length(t_target),length(t_target)]);
% score(score == 0) = NaN;
% score = score / max(max(abs(score)));
score(score > thres_noise) = 1;
score(score < -thres_noise) = -1;
score(score <= thres_noise & score >= -thres_noise) = NaN;

h = heatmap(flipud(score.'));
% h = heatmap(flipud(score));
timelabel = string(t_target);
timelabel(mod(t_target, range) ~= 0) = '';
timelabel(end) = length(t_target);
h.YDisplayLabels = flipud(timelabel);
h.XDisplayLabels = timelabel;
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';

s = struct(h);
% s.Axes.XAxisLocation = "top";
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1, 0, 1])
