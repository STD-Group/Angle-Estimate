%% ÑµÁ·¼¯
clear all;
close all;
clc;

flag = 0;

tic;

if flag
    Res = TrainGroup(2360/281600, 545/281600);
    t = toc;
    save('TrainGroup.mat');
else
    angle = TestGroup(2360/281600, 545/281600);
    t = toc;
    save('TestGroup.mat');
end

