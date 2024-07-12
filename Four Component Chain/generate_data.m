clc;
clear;
close all;

%% Parameter
R1 = random('Uniform', 3.8, 4.0);
R2 = random('Uniform', 3.5, 3.7);
R3 = random('Uniform', 3.5, 3.7);
R4 = random('Uniform', 3.7, 3.9);
A21 = random('Uniform', 0.3, 0.5);
A32 = random('Uniform', 0.3, 0.5);
A43 = random('Uniform', 0.3, 0.5);

time_step = 3000;
num_component = 4;

%% Data initialize

y = zeros(time_step+1, num_component);
t = 0:time_step;

for i = 1:num_component
    y(1,i) = random('Uniform', 0.01, 0.99);
end

for i = 1:time_step
    y(i+1, 1) = y(i, 1) * (R1 - R1 * y(i, 1));
    y(i+1, 2) = y(i, 2) * (R2 - A21 * y(i, 1) - R2 * y(i, 2));
    y(i+1, 3) = y(i, 3) * (R3 - A32 * y(i, 2) - R3 * y(i, 3));
    y(i+1, 4) = y(i, 4) * (R4 - A43 * y(i, 3) - R4 * y(i, 4));
end

%% Check for monotonicity
score = zeros(1, num_component);
score(1,1) = sum(y(:,1) > (R1-1)/(2*R1));
score(1,2) = sum(y(:,2) > (R2-1)/(2*R2));
score(1,3) = sum(y(:,3) > (R3-1)/(2*R3));
score(1,4) = sum(y(:,4) > (R4-1)/(2*R4));

save('data.mat', 'y', 't', 'score');
save('parameter.mat', 'R1', 'R2', 'R3', 'R4', 'A21', 'A32', 'A43');