
clear all;
close all;
clc;

flag = 0;

tic;

if flag
    %% ѵ����
    Res = TrainGroup(2360/281600, 545/281600);
    t = toc;
else
    %% ѵ����
    angle = TestGroup(2360/281600, 545/281600);
    t = toc;
    fid = fopen('TestGroup.txt', 'w');
    fprintf(fid, '%1.7e\n', angle);
    fclose(fid);
end

