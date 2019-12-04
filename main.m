%% ÑµÁ·¼¯
clear all;
close all;
clc;

index = 0;
k = 281600;
min = 2;
high = 0;
low = 0;
res = 0;
for High = [2350: 5: 2370]
    for Low = [520: 5: 550]
        index = index+1;
        Res(: ,  index) = TrainGroup(High/k, Low/k);
        if mean(Res(: ,  index)) < min
            min = mean(Res(: ,  index));
            high = High;
            low = Low;
        end
    end
end
res = mean(Res);