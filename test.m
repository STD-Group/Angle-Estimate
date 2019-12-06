%% ÑµÁ·¼¯
clear all;
close all;
clc;

k = 281600;

% high = [2500: 5: 2700];
high = [2200: 5: 2395];
low = [400: 5: 600];

Res = zeros(length(high), length(low), 'double');
Var = zeros(length(high), length(low), 'double');

tic;

for row = 1: length(high)
    for col = 1: length(low)
        High = high(row);
        Low = low(col);
        tmp = TrainGroup(High/k, Low/k);
        Res(row, col) = mean(abs(tmp));
        Var(row, col) = var(tmp);
    end
    save('test.mat');
end

t = toc;

save('test.mat');

