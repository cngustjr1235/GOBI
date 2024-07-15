clc;
clear;
close all;

%% load raw data
load('raw_data_repressilator')

tnorm = t / max(t);

%% plot raw data
figure(1)
plot(tnorm, y(:,1), 'Color', '#00b1f5', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
xlabel('Normalized Time')
ylabel('TetR')
yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ylim([0,1])
set(gca,'fontsize',16)


figure(2)
plot(tnorm, y(:,2), 'Color', '#ff008e', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
xlabel('Normalized Time')
ylabel('Lacl')
ylim([0,1])
yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
set(gca,'fontsize',16)


figure(3)
plot(tnorm, y(:,3), 'Color', '#ffac0a', 'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 5, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'black')
xlabel('Normalized Time')
ylabel('\lambda cl')
yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
ylim([0,1])
set(gca,'fontsize',16)
