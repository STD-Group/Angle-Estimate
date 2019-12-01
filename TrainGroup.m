clear all;
close all;
clc;
%% 文件读取
fid = fopen('train\angle.txt');
tmp = textscan(fid, '%f');
angleAns = tmp{1};
trainNum = length(angleAns);

%% 音频处理
angle = zeros(trainNum, 1, 'double');

global wave;
k = 10;
for index = 1: trainNum
    [wave, Fs] = audioread(['train\', num2str(index), '.wav']);
    wave = resample(wave, k, 1);
    Fs = Fs*k;
%     figure;
%     subplot(2, 1, 1);
%     plot(wave(: , 1));
%     subplot(2, 1, 2);
%     plot(wave(: , 2));
    angle(index) = AngleEstimate(Fs);
end

Res = sum(abs(angle-angleAns))/trainNum;
disp(['Res: ', num2str(Res)]);

%%
function angle = AngleEstimate(Fs)
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
    corrF = psdF1 .* conj(psdF2);
    SCOT = 1 ./ (psdF1 .* psdF2).^0.5;
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
    corrF11 = psdF1 .* conj(psdF1);
    corrF12 = psdF1 .* conj(psdF2);
    corrF22 = psdF2 .* conj(psdF2);
    
    corrF = corrF11 .* conj(corrF12);
    SCOT = 1 ./ (corrF11 .* corrF22).^0.5;
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
