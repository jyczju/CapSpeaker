%%本代码用于生成PWM调制波形及频谱图

clc;
clear;

time_len = 2;% 设置时间长度为2s

fs = 44100;
t_orig = 0:1/fs:time_len-1/fs;
data_read = sin(2*pi*200*t_orig); % 手动生成单频信号

filename = '200Hz_sine.mp3';

% filename = 'call-my-wife-George.mp3';
% [data_read, fs] = audioread(['voice_command_test\',filename]); % 读取音频

data_read = data_read / max(abs(data_read)); % 音频归一化

% 进行长度填充
if(length(data_read) > time_len * fs) % fs为音频采样率
    data = data_read(1:time_len * fs);
else
    data = zeros(time_len * fs, 1);
    data(1:length(data_read)) = data_read;
end

t_orig = 0:1/fs:time_len-1/fs;

figure()
plot(t_orig, data)
title("原始音频")


%% 绘制原始音频频谱图

fft_data = DrawFFT(data, fs, '原始音频频谱图');

%% AM调制
% fc = 3000;
% % data_am = (data+1)/2;
% data_am = modulate(data, fc, fs, 'am');
% 
% % figure()
% % plot(t_orig, data_am)
% % title("AM调制声音信号")
% 
% fft_data = DrawFFT(data_am, fs, 'AM调制频谱图');


%% Remove the part of original audio greater than 2000 Hz
data = mean(data, 2); % 适用于多声道，对多个声道取均值

N = length(data);
t_canon = 0 : 1/fs : (N-1)/fs; % 原始声音信号采样对应的时间点
fx_canon = 0 : fs/N : fs - fs/N;

% filter low frequency bands
pass_low = 20;
pass_high = 2000;   
ts = timeseries(data, t_canon); % 创建一个timevals对象

ts_filtered = idealfilter(ts, [pass_low, pass_high], 'pass');
data_filtered = ts_filtered.Data;

% figure()
% plot(data_filtered)
% title("低通滤波后的音频")

%% Set PWM parameters
target_frequency = 32000; % PWM carrier frequency
duty_upper_bound = 0.99; % maximum duty cycle
duty_lower_bound = 0.01; % minimum duty cycle
full_busy = 2047; % Timer accuracy decreased by 1
full_idle = 2047;% Timer accuracy decreased by 1
sample_rate = target_frequency * 100; % 对PWM波的采样率 % used to calculate the sampling rate for computing PWM duty cycle

period = 1 / target_frequency; % PWM duty period
N = time_len * target_frequency; % 对duty的采样点数
t = 0 : 1 / target_frequency : (N - 1) / target_frequency; % 对PWM波的采样时间点

%% 对原始声音信号的采样点进行变换，使其与target_frequency对齐
data_interp = interp1(t_canon', data_filtered, t', 'nearest');
target_wave = data_interp;
target_wave(find(isnan(target_wave))) = 0;

target_wave = target_wave / max(abs(target_wave)); % 再次归一化

figure()
plot(target_wave)

%% Calculate PWM duty cycle
% % 对target_wave进行非线性化
% for i = 1:length(target_wave)
%     if target_wave(i) > 0.005
%         target_wave(i) = sqrt(target_wave(i));
%     elseif target_wave(i) < -0.005
%         target_wave(i) = - sqrt(-target_wave(i));
%     end
% end
% 
% figure()
% plot(target_wave)
% title("非线性化后的声音波形")
% audiowrite(['voice_command_nonlinear/',filename(1:end-4),'_nonlinear.wav'],data,fs);
% 
% 


% 生成duty
duty = (target_wave + 1) / 2 * (duty_upper_bound - duty_lower_bound) + duty_lower_bound;

figure()
plot(duty)
title("占空比")

% fft_data = DrawFFT(duty, target_frequency, 'duty波频谱图');

busy_time = duty * full_busy; % us
results = [round(busy_time)];


%% generate pwm wave
period_pwm = 1 / sample_rate;
N_pwm = time_len * sample_rate;
t_pwm = 0 : period_pwm : (N_pwm - 1) * period_pwm;
pwm_wave = zeros(N_pwm, 1);

for i = 1:100:N_pwm-100+1
    busy_num = round(duty((i-1)/100+1)*100.0);
    pwm_wave(i:i+busy_num) = 1;  
end

figure()
plot(pwm_wave)
ylim([-1,2])
title("pwm wave")

%% 绘制PWM波频谱图
fft_data = DrawFFT(pwm_wave, sample_rate, 'PWM波频谱图');

%% Write duty cycle traces to txt files
fid = fopen(['traces_test\',filename(1:end-4),'_duty_cycle_32k.txt'], 'w');
% fid = fopen(['traces_test\',filename(1:end-4),'_duty_cycle_32k_nonlinear.txt'], 'w');
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

%% 画出信号的频谱
% data 需要处理的原始信号 fs:采样频率
function fft_data = DrawFFT(data, fs, til)
    N=length(data);
    fft_data=fft(data);
    magY=abs(fft_data(1:N/2))*2/N;
    f=(0:N/2-1)'*fs/N;
    figure()
%     plot(f(2:end),magY(2:end),'LineWidth',1.2);% 忽略直流分量
    plot(f,magY,'LineWidth',1.2);
%     h=stem(f, magY, 'fill','--');
%     set(h,'MarkerEdgeColor','red','Marker','*');
%     grid on;
    title(til);
    xlabel('f(Hz)'), ylabel('幅值');
    xlim([1,5000]) % 忽略直流分量
end