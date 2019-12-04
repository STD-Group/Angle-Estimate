%% ÑµÁ·¼¯
clear all;
close all;
clc;

index = 0;
k = 281600;
min = 2.5;
high = 0;
low = 0;
res = 0;
for High = [2590: 2: 2610]
    for Low = [520: 1: 530]
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