function Res = TrainGroup()
%     clear all;
%     close all;
%     clc;
    %% 文件读取
    fid = fopen('train\angle.txt');
    tmp = textscan(fid, '%f');
    angleAns = tmp{1};
    trainNum = length(angleAns);

    %% 音频处理
    angle = zeros(trainNum, 1, 'double');
    
    % 过拟合得到的参数
    High = 2600 / 281600;
    Low = 524 / 281600;
    
    global wave;
    k = 10;
    for index = 1: trainNum
        [wave, Fs] = audioread(['train\', num2str(index), '.wav']);
        wave = resample(wave, k, 1);
        Fs = Fs*k;

        angle(index) = AngleEstimate(Fs, High, Low);
    end
    Res = abs(angle-angleAns);
    disp(['Res: ', num2str(sum(Res)/trainNum)]);
    disp(['Variance: ', num2str(var(Res))]);
end

%%
function angle = AngleEstimate(Fs, High, Low)
    speed = 343;
    global wave;
%     psd = pwelch(wave);
%     SCOT = 1 / sqrt(conv(psd(: , 1), psd(: , 2)));
%     corrF = xcorr(psd(: , 1), psd(: , 2));
%     [~, M] = max(corrF);
%     disp(M / length(corrF) - 0.5);
    
    L = size(wave, 1);
    psdF = fft(wave);
    psdF1 = psdF(: , 1);
    psdF2 = psdF(: , 2);
    
    %% 对功率谱FFT结果进行滤波
    FilterF = ones(length(psdF1), 1, 'double')*0.001;
    indexHigh = round(High*length(FilterF));
    indexLow = round(Low*length(FilterF));
%     tmp = abs(psdF1(indexLow: indexHigh));
%     [~, peak] = max(tmp);
%     peak = peak+indexLow-1;
%     indexHigh = round(peak * alphaHigh);
%     indexLow = round(peak * alphaLow);
    FilterF(indexLow: indexHigh) = 1;
    FilterF(end-indexHigh+1: end-indexLow+1) = 1;
    psdF1 = psdF1 .* FilterF;
    psdF2 = psdF2 .* FilterF;
    
    
    %% 计算相关函数
    corrF = psdF1 .* conj(psdF2);
    SCOT = 1 ./ abs((psdF1 .* psdF2).^0.5);
    PHAT = 1 ./ abs(corrF);
%     corrT = ifft(corrF .* PHAT);
%     corrT = ifft(corrF .* SCOT);
    corrT = ifft(corrF);
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
    
    disp(['angle: ', num2str(angle)]);

end

function angle = SecondaryCorrelation(Fs)
    speed = 343;
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
    
%     corrT = ifft(corrF);
%     corrT = ifft(corrF .* PHAT);
    corrT = ifft(corrF .* SCOT);
%     plot(corrT)

%     SCOT = 1 ./ abs(psdF1 .* conj(psdF2)).^0.65;
%     PHAT = 1 ./ abs(corrF12);
%     corrT = abs(ifft(corrF .* PHAT));
%     corrT = abs(ifft(corrF .* SCOT));

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
    disp(['angle: ', num2str(angle)]);

end
