function Res = TrainGroup(High, Low)
    %% 文件读取
    fid = fopen('train\angle.txt');
    tmp = textscan(fid, '%f');
    angleAns = tmp{1};
    trainNum = length(angleAns);
    fclose(fid);

    %% 音频处理
    angle = zeros(trainNum, 1, 'double');
    
    global wave;
    k = 10;
    for index = 1: trainNum
        [wave, Fs] = audioread(['train\', num2str(index), '.wav']);
        wave = resample(wave, k, 1);
        Fs = Fs*k;
        % 主要处理函数
        angle(index) = AngleEstimate(Fs, High, Low);
    end
    Res = abs(angle-angleAns);
%     disp(['Res: ', num2str(sum(Res)/trainNum)]);
%     disp(['Variance: ', num2str(var(Res))]);

end

%%
function angle = AngleEstimate(Fs, High, Low)
    speed = 343;
    strategy = 0; % 窗函数选择
    global wave;
    
    L = size(wave, 1);
    psdF = fft(wave);
    psdF1 = psdF(: , 1);
    psdF2 = psdF(: , 2);
    
    %% 对功率谱FFT结果进行滤波
    FilterF = ones(length(psdF1), 1, 'double')*0.001;
    indexHigh = round(High*length(FilterF));
    indexLow = round(Low*length(FilterF));
    FilterF(indexLow: indexHigh) = 1;
    FilterF(end-indexHigh+1: end-indexLow+1) = 1;
    psdF1 = psdF1 .* FilterF;
    psdF2 = psdF2 .* FilterF;
    
    %% 计算相关函数
    corrF = psdF1 .* conj(psdF2);
    SCOT = 1 ./ abs((psdF1 .* psdF2).^0.5);
    PHAT = 1 ./ abs(corrF);
    switch strategy
        case 1
            corrT = ifft(corrF .* SCOT);
        case 2
            corrT = ifft(corrF .* PHAT);
        otherwise
            corrT = ifft(corrF);
    end

%     plot(corrT)
    %% 确定时域相关峰值
    len = 0.1;
    [~, M] = max(corrT);
    if M <= L / 2
        dis = M / Fs * speed;
        if dis <= len
            angle = acos(dis/len)*180 / pi;
        else
            angle = 0;
        end
    else
        dis = (L-M) / Fs * speed;
        if dis <= len
            angle = 180 -acos(dis/len)*180 / pi;
        else
            angle = 180;
        end
    end
    
%     disp(['angle: ', num2str(angle)]);

end

function angle = SecondaryCorrelation(Fs)
    speed = 343;
    strategy = 0; % 窗函数选择
    global wave;
    
    L = size(wave, 1);
    psdF = fft(wave);
    psdF1 = psdF(: , 1);
    psdF2 = psdF(: , 2);
%% 对功率谱FFT结果进行滤波
    FilterF = ones(length(psdF1), 1, 'double')*0.01;;
    indexHigh = 2500;
    indexLow = 200;
    FilterF(indexLow: indexHigh) = 1;
    FilterF(end-indexHigh: end-indexLow) = 1;
    psdF1 = psdF1 .* FilterF;
    psdF2 = psdF2 .* FilterF;
    
%% 计算相关函数
    corrF11 = psdF1 .* conj(psdF1);
    corrF12 = psdF1 .* conj(psdF2);
    corrF22 = psdF2 .* conj(psdF2);
    
    corrF = corrF11 .* conj(corrF12);
    SCOT = 1 ./ abs((corrF11 .* corrF22).^0.5);
    PHAT = 1 ./ abs(corrF);
    
    switch strategy
        case 1
            corrT = ifft(corrF .* SCOT);
        case 2
            corrT = ifft(corrF .* PHAT);
        otherwise
            corrT = ifft(corrF);
    end
    
%     plot(corrT)

    len = 0.1;
    [~, M] = max(corrT);
    if M <= L / 2
        dis = M / Fs * speed;
        if dis <= len
            angle = acos(dis/len)*180 / pi;
        else
            angle = 0;
        end
    else
        dis = (L-M) / Fs * speed;
        if dis <= len
            angle = 180 -acos(dis/len)*180 / pi;
        else
            angle = 180;
        end
    end
    
    angle = 180-angle;
%     disp(['angle: ', num2str(angle)]);

end
