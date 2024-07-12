clc;
clear;
close all;
addpath('../');
addpath('GOBI');
addpath('../Utility');

%% load raw data
load('data_ori.mat')
t_ori = data_ori.DayNumber;
y_ori = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
labels = ["Cyclopoids", "Protozoa", "Rotifers"];

num_component = length(y_ori(1,:));

y_ori(isnan(y_ori)) = 0;

%% options 1. choose thresholds and critical values

% defaults
noise_level = 30; % compute noise level using residuals
thres_R = 0;     % threshold for regulation-detection region
                 % thres_R = 0 represents default (decreases as dimension increases)
                 
thres_S = 0.9 - 0.005 * noise_level; % threshold for regulation-detection score
% thres_S = 0.2;
thres_TRS = 0.9 - 0.01 * noise_level;% threshold for total regulation score

p_delta = 0.01; % critical value for delta test
p_surrogate = 0.01; % critical value for surrogate test

%% options 2. Choose the types of self-regulations

% type_self = nan: infer without any assumptions of self regulations
% type_self = -1 : negative self regulations
% type_self = 0  : no self regulations
% type_self = 1  : positive self regulations

type_self = nan;

%% options 3. maximum dimension of framework

% default
if isnan(type_self)
    max_D = num_component;
else
    max_D = num_component - 1;
end

%% equidistant
spline_param = 0.01;

t_eq = (t_ori(1):1:t_ori(end)).';
y_eq = zeros(length(t_eq), num_component);
for i = 1:num_component
    y_eq(:, i) = csaps(t_ori, y_ori(:, i), spline_param, t_eq);
end

%% run with delay - total series

y_delay_total = {};
y_fit_total = {};
t_fit_total = {};
residual_total = {};

S_total_list_1D_total = {};
R_total_list_1D_total = {};
TRS_list_1D_total = {};
detection_list_1D_total = {};

S_total_list_2D_total = {};
R_total_list_2D_total = {};
TRS_list_2D_total = {};
detection_list_2D_total = {};

delta_total = {};
surrogate_total = {};
regulation_result_total = {};
fixed_target = 1;

for delay_pro = 0:10
    for delay_rot = 0:10
        % Preprocess data - cut and interpolate

        max_delay = max(delay_pro, delay_rot);
        
        cyc = y_eq(1+max_delay:end,1);
        pro = y_eq(1+max_delay-delay_pro:end-delay_pro,2);
        rot = y_eq(1+max_delay-delay_rot:end-delay_rot,3);

        y_delay = [cyc, pro, rot];

        length_timeseries = length(cyc);
        window_size = length_timeseries;
        time_interval = 1/(window_size-1);

        [y_fit, t_fit, num_data, residual] = cut_interpolate(y_delay, length_timeseries, 0, "none");

        y_delay_total{delay_pro+1, delay_rot+1} = y_delay;
        y_fit_total{delay_pro+1, delay_rot+1} = y_fit;
        t_fit_total{delay_pro+1, delay_rot+1} = t_fit;
        residual_total{delay_pro+1, delay_rot+1} = residual;

        % Compute 1D RDS
        thres_noise_1D = 0;
        [S_total_list_1D, R_total_list_1D, pair_list_1D] = RDS_1D(y_fit, t_fit, fixed_target, num_component, type_self, time_interval, thres_noise_1D);
        S_total_list_1D_total{delay_pro+1, delay_rot+1} = S_total_list_1D;
        R_total_list_1D_total{delay_pro+1, delay_rot+1} = R_total_list_1D;

        % Compute 1D TRS
        if thres_R == 0
            thres_R = 0.05;
        end
        [TRS_list_1D, detection_list_1D] = TRS(S_total_list_1D, R_total_list_1D, pair_list_1D, thres_S, thres_R, thres_TRS, 1);
        TRS_list_1D_total{delay_pro+1, delay_rot+1} = TRS_list_1D;
        detection_list_1D_total{delay_pro+1, delay_rot+1} = detection_list_1D;

        % Compute 2D RDS
        thres_noise_2D = 0;
        [S_total_list_2D, R_total_list_2D, pair_list_2D] = RDS_2D(y_fit, t_fit, fixed_target, num_component, type_self, time_interval, thres_noise_2D);
        S_total_list_2D_total{delay_pro+1, delay_rot+1} = S_total_list_2D;
        R_total_list_2D_total{delay_pro+1, delay_rot+1} = R_total_list_2D;

        % Compute 2D TRS
        if thres_R == 0
            thres_R = 0.05;
        end
        [TRS_list_2D, detection_list_2D] = TRS(S_total_list_2D, R_total_list_2D, pair_list_2D, thres_S, thres_R, thres_TRS, 2);
        TRS_list_2D_total{delay_pro+1, delay_rot+1} = TRS_list_2D;
        detection_list_2D_total{delay_pro+1, delay_rot+1} = detection_list_2D;

        % Delta Test
        [delta_candidate_list_2D, delta_type_list_2D, delta_list_2D] = delta_test_2D(S_total_list_2D, R_total_list_2D, detection_list_2D, pair_list_2D, thres_R);
        delta_total{delay_pro+1, delay_rot+1} = [delta_candidate_list_2D, delta_type_list_2D, delta_list_2D];

        % Surrogate Test
        [boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D] = surrogate_test_2D(y_fit, t_fit, delta_candidate_list_2D, delta_type_list_2D, delta_list_2D, p_delta, p_surrogate, type_self, thres_noise_2D);
        surrogate_total{delay_pro+1, delay_rot+1} = [boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D];

        % Merge regulation
        [regulation_network_1D, regulation_network_2D] = merge_regulation(pair_list_1D, detection_list_1D, pair_list_2D, boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D, num_component);
        regulation_result_total{delay_pro+1, delay_rot+1} = regulation_network_2D;

    end
end

save('Delay_cut');

%% run with delay - cut series

y_delay_total = {};
y_fit_total = {};
t_fit_total = {};
residual_total = {};

S_total_list_1D_total = {};
R_total_list_1D_total = {};
TRS_list_1D_total = {};
detection_list_1D_total = {};

S_total_list_2D_total = {};
R_total_list_2D_total = {};
TRS_list_2D_total = {};
detection_list_2D_total = {};

delta_total = {};
surrogate_total = {};
regulation_result_total = {};
fixed_target = 1;

for delay_pro = 0:10
    for delay_rot = 0:10
        % Preprocess data - cut and interpolate

        max_delay = max(delay_pro, delay_rot);
        
        y1 = y_eq(216:792, :);
        cyc1 = y1(1+max_delay:end,1);
        pro1 = y1(1+max_delay-delay_pro:end-delay_pro,2);
        rot1 = y1(1+max_delay-delay_rot:end-delay_rot,3);
        y_delay1 = [cyc1, pro1, rot1];

        y2 = y_eq(1027:end, :);
        cyc2 = y2(1+max_delay:end,1);
        pro2 = y2(1+max_delay-delay_pro:end-delay_pro,2);
        rot2 = y2(1+max_delay-delay_rot:end-delay_rot,3);
        y_delay2 = [cyc2, pro2, rot2];

        length_timeseries = length(y_delay1(:,1));
        window_size = length_timeseries;
        time_interval = 1/(window_size-1);

        [y_fit1, t_fit1, num_data1, residual1] = cut_interpolate(y_delay1, length_timeseries, 0.9);
        % cut second part with window size as first part length
        % t_fit1 is equal to t_fit2
        [y_fit2, t_fit2, num_data2, residual2] = cut_interpolate(y_delay2, length_timeseries, 0.9);

        y_fit = [y_fit1, y_fit2];
        num_data = num_data1 + num_data2;
        t_fit = t_fit1;

        y_delay_total{delay_pro+1, delay_rot+1} = {y_delay1, y_delay2};
        y_fit_total{delay_pro+1, delay_rot+1} = y_fit;
        t_fit_total{delay_pro+1, delay_rot+1} = t_fit;
        % residual_total{delay_pro+1, delay_rot+1} = residual;

        % Compute 1D RDS
        thres_noise_1D = 0;
        [S_total_list_1D, R_total_list_1D, pair_list_1D] = RDS_1D(y_fit, t_fit, fixed_target, num_component, type_self, time_interval, thres_noise_1D);
        S_total_list_1D_total{delay_pro+1, delay_rot+1} = S_total_list_1D;
        R_total_list_1D_total{delay_pro+1, delay_rot+1} = R_total_list_1D;

        % Compute 1D TRS
        if thres_R == 0
            thres_R = 0.05;
        end
        [TRS_list_1D, detection_list_1D] = TRS(S_total_list_1D, R_total_list_1D, pair_list_1D, thres_S, thres_R, thres_TRS);
        TRS_list_1D_total{delay_pro+1, delay_rot+1} = TRS_list_1D;
        detection_list_1D_total{delay_pro+1, delay_rot+1} = detection_list_1D;

        % Compute 2D RDS
        thres_noise_2D = 0;
        [S_total_list_2D, R_total_list_2D, pair_list_2D] = RDS_2D(y_fit, t_fit, fixed_target, num_component, type_self, time_interval, thres_noise_2D);
        S_total_list_2D_total{delay_pro+1, delay_rot+1} = S_total_list_2D;
        R_total_list_2D_total{delay_pro+1, delay_rot+1} = R_total_list_2D;

        % Compute 2D TRS
        if thres_R == 0
            thres_R = 0.05;
        end
        [TRS_list_2D, detection_list_2D] = TRS(S_total_list_2D, R_total_list_2D, pair_list_2D, thres_S, thres_R, thres_TRS);
        TRS_list_2D_total{delay_pro+1, delay_rot+1} = TRS_list_2D;
        detection_list_2D_total{delay_pro+1, delay_rot+1} = detection_list_2D;

        % Delta Test
        [delta_candidate_list_2D, delta_type_list_2D, delta_list_2D] = delta_test_2D(S_total_list_2D, R_total_list_2D, detection_list_2D, pair_list_2D, thres_R);
        delta_total{delay_pro+1, delay_rot+1} = [delta_candidate_list_2D, delta_type_list_2D, delta_list_2D];

        % Surrogate Test
        [boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D] = surrogate_test_2D(y_fit, t_fit, delta_candidate_list_2D, delta_type_list_2D, delta_list_2D, p_delta, p_surrogate, type_self, thres_noise_2D);
        surrogate_total{delay_pro+1, delay_rot+1} = [boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D];

        % Merge regulation
        [regulation_network_1D, regulation_network_2D] = merge_regulation(pair_list_1D, detection_list_1D, pair_list_2D, boot_candidate_list_2D, boot_type_list_2D, surrogate_list_2D, num_component);
        regulation_result_total{delay_pro+1, delay_rot+1} = regulation_network_2D;

    end
end

save('Delay_cut_120_total');


