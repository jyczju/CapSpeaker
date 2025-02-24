%% Read the original audio, and the audio length is complemented to 2 seconds.
clc;clear;
tic;

time_len = 2;
[data_read, fs] = audioread('voice_command_test\open-the-door-George.mp3'); % 读取音频

% figure()
% plot(data_read)

data_read = data_read / max(abs(data_read)); % 音频归一化

% 进行长度填充
if(length(data_read) > time_len * fs) % fs为采样率
    data = data_read(1:time_len * fs);
else
    data = zeros(time_len * fs, 1);
    data(1:length(data_read)) = data_read;
end

% figure()
% plot(data)

%% Remove the part of original audio greater than 2000 Hz
data = mean(data, 2); % 适用于多声道，对多个声道取均值

N = length(data);
t_canon = 0 : 1/fs : (N-1)/fs; % 采样对应的时间点
fx_canon = 0 : fs/N : fs - fs/N; % 不知道干嘛的
data_fft = abs(fft(data)); % 对时间数据作傅里叶变换

% figure()
% plot(data_fft)

% filter low frequency bands
pass_low = 20;
pass_high = 2000;   
ts = timeseries(data, t_canon); % 创建一个timevals对象

ts_filtered = idealfilter(ts, [pass_low, pass_high], 'pass');
data_filtered = ts_filtered.Data;

% figure()
% plot(data_filtered)

%% Set PWM parameters
target_frequency = 24000; % PWM carrier frequency
duty_upper_bound = 0.99; % maximum duty cycle
duty_lower_bound = 0.01; % minimum duty cycle
full_busy = 2047; % Timer accuracy decreased by 1
full_idle = 2047;% Timer accuracy decreased by 1
sample_rate = target_frequency * 5; % used to calculate the sampling rate for computing PWM duty cycle

period = 1 / target_frequency; % PWM period
N = time_len * sample_rate;
t = 0 : 1 / sample_rate : (N - 1) / sample_rate;

%% The original audio is interpolated for easy calculation
data_interp = interp1(t_canon', data_filtered, t', 'nearest');
target_wave = data_interp;
target_wave(find(isnan(target_wave))) = 0;

% figure()
% plot(target_wave)

%% Calculate PWM duty cycle
period_samples = round(period * sample_rate); % =5
period_num = floor(N / period_samples);
duty = zeros(period_num, 1);
pwm_wave = zeros(N, 1);
for ii = 1:period_num
    strength = target_wave((ii-1) * period_samples + 1 ); % 隔5个取1个
    duty(ii) = (strength + 1) / 2 * (duty_upper_bound - duty_lower_bound) + duty_lower_bound;
    busy_samples = round(duty(ii) * period_samples);
    pwm_wave((ii-1) * period_samples + 1: (ii-1) * period_samples + busy_samples) = 1;
end
busy_time = duty * full_busy; % us
idle_time = (1 - duty) * full_idle;
results = [round(busy_time)];
%% Write duty cycle traces to txt files
fid = fopen(['traces_test\copen-the-door-George_duty_cycle_24k.txt'], 'w');
fprintf(fid, "a={");
for i=1:2:length(results)
    fprintf(fid, "%d,", results(i));   
    if(mod(i, 20)== 1) 
        fprintf(fid, "\n");
    end
end
fprintf(fid, "%d", results(end));
fprintf(fid, "};");
fclose(fid);
    

toc