%% Protozoa Dominant Phase

load('data_ori');

t = data_ori.DayNumber;
y = [data_ori.Cyclopoids, data_ori.Protozoa, data_ori.Rotifers];
t = t(79:end);
y = y(79:end, :);

%% Phase 1

t1 = t(66:88);
y1 = y(66:88, :);

y1(:, 1) = y1(:, 1) - min(y1(:, 1));
if max(y1(:, 1)) ~= 0
    y1(:, 1) = y1(:, 1) / max(y1(:, 1));
end

figure(Position=[0, 0, 400, 300]);
colororder({'k', 'k'});
fontsize(20, 'points');
yyaxis left;
% plot(t1, y1(:,1), Color="red", LineWidth=2, LineStyle="-");
plot(t1, y1(:,2), Color="#2A8947", LineWidth=2, LineStyle="-");
hold on;
yyaxis right
plot(t1, y1(:,3), Color="blue", LineWidth=2, LineStyle="-");
xlim([t1(1), t1(end)]);
legend(["Protozoa", "Rotifers"]);

yyaxis right
ylabel("Rotifers")
yyaxis left
ylabel("Protozoa")

%% Phase 2

t2 = t(322:387);
y2 = y(322:387, :);

y2(:, 1) = y2(:, 1) - min(y2(:, 1));
if max(y2(:, 1)) ~= 0
    y2(:, 1) = y2(:, 1) / max(y2(:, 1));
end

figure(Position=[0, 0, 400, 300]);
colororder({'k', 'k'});
yyaxis left;

plot(t2, y2(:,2), Color="#2A8947", LineWidth=2, LineStyle="-");
hold on;
yyaxis right
plot(t2, y2(:,3), Color="blue", LineWidth=2, LineStyle="-");
xlim([t2(1), t2(end)]);
legend(["Protozoa", "Rotifers"]);

yyaxis right
ylabel("Rotifers")
yyaxis left
ylabel("Protozoa")

%% Phase 3

t3 = t(576:625);
y3 = y(576:625, :);

y3(:, 1) = y3(:, 1) - min(y3(:, 1));
if max(y3(:, 1)) ~= 0
    y3(:, 1) = y3(:, 1) / max(y3(:, 1));
end

figure(Position=[0, 0, 400, 300]);
colororder({'k', 'k'});
yyaxis left;
plot(t3, y3(:,2), Color="#2A8947", LineWidth=2, LineStyle="-");
hold on;
yyaxis right
plot(t3, y3(:,3), Color="blue", LineWidth=2, LineStyle="-");
xlim([t3(1), t3(end)]);
legend(["Protozoa", "Rotifers"]);

yyaxis right
ylabel("Rotifers")
yyaxis left
ylabel("Protozoa")

%% Rotifers Dominant Phase
%% Phase 1

t1 = t(71:106);
y1 = y(71:106, :);

y1(:, 1) = y1(:, 1) - min(y1(:, 1));
if max(y1(:, 1)) ~= 0
    y1(:, 1) = y1(:, 1) / max(y1(:, 1));
end

figure(1);
colororder({'k', 'k'});
yyaxis left;
plot(t1, y1(:,1), Color="red", LineWidth=2, LineStyle="-");
hold on;
plot(t1, y1(:,2), Color="#2A8947", LineWidth=2, LineStyle="-");
yyaxis right
plot(t1, y1(:,3), Color="blue", LineWidth=2, LineStyle="-");
xlim([t1(1), t1(end)]);
legend(["Normalized Cyclopoids", "Protozoa", "Rotifers"]);

yyaxis right
ylabel("Rotifers")
yyaxis left
ylabel("Cyclopoids and Protozoa")

%% Phase 2

t2 = t(636:685);
y2 = y(636:685, :);

y2(:, 1) = y2(:, 1) - min(y2(:, 1));
if max(y2(:, 1)) ~= 0
    y2(:, 1) = y2(:, 1) / max(y2(:, 1));
end

figure(2);
colororder({'k', 'k'});
yyaxis left;
plot(t2, y2(:,1), Color="red", LineWidth=2, LineStyle="-");
hold on;
plot(t2, y2(:,2), Color="#2A8947", LineWidth=2, LineStyle="-");
yyaxis right
plot(t2, y2(:,3), Color="blue", LineWidth=2, LineStyle="-");
xlim([t2(1), t2(end)]);
legend(["Normalized Cyclopoids", "Protozoa", "Rotifers"]);

yyaxis right
ylabel("Rotifers")
yyaxis left
ylabel("Cyclopoids and Protozoa")
