clc;
clear;

time_len = 1.0;% 设置时间长度为2s

filename = '31k-sqrt-iPhone4s.wav';
[data_read, fs] = audioread(['voice_plot\',filename]); % 读取音频

% 进行长度填充
if(length(data_read) > time_len * fs) % fs为音频采样率
    data = data_read(1:time_len * fs);
else
    data = zeros(time_len * fs, 1);
    data(1:length(data_read)) = data_read;
end

t = 0:1/fs:time_len-1/fs;

figure()
subplot(1,2,1)
plot(t, data)
ylim([-0.003,0.003])
title("Calibration Audio")

filename = '31k-linear-iPhone4s.wav';
[data_read, fs] = audioread(['voice_plot\',filename]); % 读取音频

% 进行长度填充
if(length(data_read) > time_len * fs) % fs为音频采样率
    data = data_read(1:time_len * fs);
else
    data = zeros(time_len * fs, 1);
    data(1:length(data_read)) = data_read;
end

subplot(1,2,2)
plot(t, data)
ylim([-0.003,0.003])
title("Calibration Audio")