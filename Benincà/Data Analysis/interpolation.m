clc;
clear;
close all;

%% load raw data
load('data.mat')
t = data.DayNumber;
y = horzcat(data.Rotifers, data.Calanoids, data.Picophytoplankton, data.Nanophytoplankton);
label = ["Rotifers", "Calanoids", "Pico", "Nano"];
color = ["#D95319", "#7E2F8E", "#0072BD", "#77AC30"];


%% parameters

num_component = length(y(1,:));
window_size_ori = 90;
time_interval = 1/10;

t_fit = linspace(0, t(end), length(t)/time_interval + 1).';

%% overall plot

figure("Position",[0, 0, 1000, 500])
for i = 1:num_component
    plot(t, y(:,i), Color=color(i))
    hold on
end
xlim([0, t(end)])
ylim([0, max(y, [], 'all')])
legend(label(:))
title('Overall', FontSize=16)

%% interpolate - spline
y_spline = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_spline(:, i) = interp1(t, y(:, i), t_fit, 'spline');
end

%% spline plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])
    plot(t, y(:, i),'k')
    hold on
    plot(t_fit, y_spline(:, i),Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y(:, i)), max(y_spline(:, i)))])
    legend('data', 'spline')
    title(['Splined ', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y(:, i),'k-o')
    hold on
    plot(t_fit, y_spline(:, i), Color=color(i),Marker='o');
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y(1:window_size_ori, i)), max(y_spline(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'spline')
    title(['Splined in Window ', label(i)], FontSize=16)
end

%% interpolate - fourier4
y_fourier4 = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_fourier4(:, i) = interp1(t, y(:, i), t_fit, 'linear');
    fitting = fit(t_fit, y_fourier4(:,i), 'fourier4');
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for j = 1:4
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_fourier4(:, i) = coeffs(1:end-1) * fouriers;
end

%% fourier plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])   
    plot(t, y(:, i),'k')
    hold on
    plot(t_fit, y_fourier4(:, i), Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y(:, i)), max(y_fourier4(:, i)))])
    legend('data', 'fourier4')
    title(['Fourier4', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y(:, i), 'k-o')
    hold on
    plot(t_fit, y_fourier4(:, i), Color=color(i), Marker='o')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y(1:window_size_ori, i)), max(y_fourier4(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'fourier4')
    title(['Fourier4 in Window', label(i)], FontSize=16)
end

%% compressed by 1/4

y_comp = nthroot(y, 4);
figure("Position",[0, 0, 1000, 500])   
for i = 1:num_component
    plot(t, y_comp(:,i), Color=color(i))
    hold on
end
xlim([0, t(end)])
ylim([0, max(y_comp, [], 'all')])
legend(label(:))
title('Compressed by 1/4', FontSize=16)

%% compressed interpolate - spline
y_comp_spline = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_comp_spline(:, i) = interp1(t, y_comp(:, i), t_fit, 'spline');
end

%% compressed spline plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])   
    plot(t, y_comp(:, i),'k')
    hold on
    plot(t_fit, y_comp_spline(:, i), Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_spline(:, i)))])
    legend('data', 'spline')
    title(['Compressed Splined', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y_comp(:, i), 'k-o')
    hold on
    plot(t_fit, y_comp_spline(:, i), Color=color(i), Marker='o')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y_comp(1:window_size_ori, i)), max(y_comp_spline(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'spline')
    title(['Compressed Splined in Window', label(i)], FontSize=16)
end

%% compressed interpolate - fourier4
y_comp_fourier4 = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_comp_fourier4(:, i) = interp1(t, y_comp(:, i), t_fit, 'linear');
    fitting = fit(t_fit, y_comp_fourier4(:,i), 'fourier4');
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for j = 1:4
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_comp_fourier4(:, i) = coeffs(1:end-1) * fouriers;
end

%% compressed fourier plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])  
    plot(t, y_comp(:, i),'k')
    hold on
    plot(t_fit, y_comp_fourier4(:, i),Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_fourier4(:, i)))])
    legend('data', 'fourier4')
    title(['Compressed Fourier4', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y_comp(:, i), 'k-o')
    hold on
    plot(t_fit, y_comp_fourier4(:, i), Color=color(i), Marker='o')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y_comp(1:window_size_ori, i)), max(y_comp_fourier4(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'fourier4')
    title(['Compressed Fourier4 in Window', label(i)], FontSize=16)
end

%% smoothing

y_smth = zeros(length(t), num_component);

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])
    y_smth(:,i) = smooth(y(:, i));
    plot(t, y(:,i),Color='k')
    hold on
    plot(t, y_smth(:,i), Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y(:, i)), max(y_smth(:, i)))])
    legend('data', 'smooth')
    title(['Smoothing', label(i)], FontSize=16)
end

%% compressed smoothing

y_comp_smth = zeros(length(t), num_component);

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])
    y_comp_smth(:,i) = smooth(y_comp(:, i));
    plot(t, y_comp(:,i),Color='k')
    hold on
    plot(t, y_comp_smth(:,i), Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_smth(:, i)))])
    legend('data', 'smooth')
    title(['Compressed Smoothing', label(i)], FontSize=16)
end

%% compressed smoothing interpolate - fourier4
y_comp_smth_fourier4 = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_comp_smth_fourier4(:, i) = interp1(t, y_comp_smth(:, i), t_fit, 'linear');
    fitting = fit(t_fit, y_comp_smth_fourier4(:,i), 'fourier4');
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for j = 1:4
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_comp_smth_fourier4(:, i) = coeffs(1:end-1) * fouriers;
end

%% compressed smoothing fourier plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])  
    plot(t, y_comp(:, i),'k')
    hold on
    plot(t_fit, y_comp_fourier4(:, i), 'k')
    hold on
    plot(t, y_comp_smth(:, i), Color=color(i))
    hold on
    plot(t_fit, y_comp_smth_fourier4(:, i),Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_smth_fourier4(:, i)))])
    legend('compressed data', 'compressed fourier4', 'compressed smoothing', 'compressed smoothing fourier4')
    title(['Compressed Smoothing Fourier4', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y_comp(:, i), 'k-o')
    hold on
    plot(t_fit, y_comp_fourier4(:, i), 'k-*')
    hold on
    plot(t, y_comp_smth(:, i), Color=color(i), Marker='o')
    hold on
    plot(t_fit, y_comp_smth_fourier4(:, i), Color=color(i), Marker='*')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y_comp(1:window_size_ori, i)),max(y_comp_smth_fourier4(1:(window_size_ori/time_interval + 1), i)))])
    legend('compressed data', 'compressed fourier4', 'compressed smoothing', 'compressed smoothing fourier4')
    title(['Compressed Smoothing Fourier4 in Window', label(i)], FontSize=16)
end

%% compressed interpolate - fourier6
y_comp_fourier6 = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_comp_fourier6(:, i) = interp1(t, y_comp(:, i), t_fit, 'linear');
    fitting = fit(t_fit, y_comp_fourier6(:,i), 'fourier6');
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for j = 1:6
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_comp_fourier6(:, i) = coeffs(1:end-1) * fouriers;
end

%% compressed fourier6 plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])  
    plot(t, y_comp(:, i),'k')
    hold on
    plot(t_fit, y_comp_fourier6(:, i),Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_fourier6(:, i)))])
    legend('data', 'fourier6')
    title(['Compressed Fourier6', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y_comp(:, i), 'k-o')
    hold on
    plot(t_fit, y_comp_fourier6(:, i), Color=color(i), Marker='o')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y_comp(1:window_size_ori, i)), max(y_comp_fourier6(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'fourier6')
    title(['Compressed Fourier6 in Window', label(i)], FontSize=16)
end

%% compressed interpolate - fourier8
y_comp_fourier8 = zeros(length(t_fit), num_component);

for i = 1:num_component
    y_comp_fourier8(:, i) = interp1(t, y_comp(:, i), t_fit, 'linear');
    fitting = fit(t_fit, y_comp_fourier8(:,i), 'fourier8');
    w = fitting.w;
    fouriers = ones(1, length(t_fit));
    for j = 1:8
        fouriers = [fouriers; cos(j*w*t_fit.')];
        fouriers = [fouriers; sin(j*w*t_fit.')];
    end
    coeffs = coeffvalues(fitting);
    y_comp_fourier8(:, i) = coeffs(1:end-1) * fouriers;
end

%% compressed fourier8 plot

for i = 1:num_component
    figure("Position",[0, 0, 1000, 500])  
    plot(t, y_comp(:, i),'k')
    hold on
    plot(t_fit, y_comp_fourier8(:, i),Color=color(i))
    xlim([0, t(end)])
    ylim([0, max(max(y_comp(:, i)), max(y_comp_fourier8(:, i)))])
    legend('data', 'fourier8')
    title(['Compressed Fourier8', label(i)], FontSize=16)

    figure("Position",[0, 0, 1000, 500])
    plot(t, y_comp(:, i), 'k-o')
    hold on
    plot(t_fit, y_comp_fourier8(:, i), Color=color(i), Marker='o')
    xlim([0, t(window_size_ori)])
    ylim([0, max(max(y_comp(1:window_size_ori, i)), max(y_comp_fourier8(1:(window_size_ori/time_interval + 1), i)))])
    legend('data', 'fourier8')
    title(['Compressed Fourier8 in Window', label(i)], FontSize=16)
end