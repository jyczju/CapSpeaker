%% Read the original audio, and the audio length is complemented to 5 seconds.
clc;clear;
tic;
filelist = dir('voice_command') % input dir


% for ff = 3:length(filelist) %filelist前两个文件名为.和..\
for ff = 5:5
    %% Audio file reading and amplifying
    filename = filelist(ff).name;
    fprintf(filename, "%d,")
    time_len = 2; %音频时长为两秒
    [data_read, fs] = audioread(['voice_command\',filename]);
%     plot(data_read);
    data_read = data_read / max(abs(data_read)); %归一化
    if(length(data_read) > time_len * fs)  % 取设定时间
         data = data_read(1:time_len * fs);
    else
        data = zeros(time_len * fs, 1);
        data(1:length(data_read)) = data_read;
    end
%     figure(1)
%     plot(data);
%     title('Original signal');

    %% 设置AM调制信号
    fc = 3000;
    data = modulate(data, fc, fs, 'am');
    figure(2)
    plot(data);
    title('AM modulated signal')
   
      %% 看信号单侧频谱分布
    figure(3);
    N1 = time_len * fs;
    y1 = abs(fft(data)/N1);
    y2 = y1(1:N1/2+1)
    f1 = (0:round(N1/2))*fs/N1;
    plot(f1,y2)
    title('Spectrogram of AM modulated signal')
%     xlim([6 inf]);

    %% Remove the part of original audio greater than 2000 Hz
    data = mean(data, 2);
    N = length(data);
    t_canon = 0 : 1/fs : (N-1)/fs;
    data_fft = abs(fft(data));
    % filter low frequency bands
    pass_low = 20;
    pass_high = 2000;   
    ts = timeseries(data, t_canon);
    ts_filtered = idealfilter(ts, [pass_low, pass_high], 'pass');
    data_filtered = ts_filtered.Data;
%     % show the original signal and filtered signal
%     plot(ts)
%     hold on
%     plot(ts_filtered);


    %% Set PWM parameters
    target_frequency = 32000; % PWM carrier frequency
    duty_upper_bound = 0.99; % maximum duty cycle
    duty_lower_bound = 0.01; % minimum duty cycle
    full_busy = 2047; % Timer accuracy decreased by 1
    full_idle = 2047;% Timer accuracy decreased by 1
    sample_rate = target_frequency *2; % used to calculate the sampling rate for computing PWM duty cycle

    period = 1 / target_frequency; % PWM period
    N = time_len * sample_rate;
    t = 0 : 1 / sample_rate : (N - 1) / sample_rate;

    %% The original audio is interpolated for easy calculation
    data_interp = interp1(t_canon', data_filtered, t', 'nearest');
    target_wave = data_interp;
    target_wave(find(isnan(target_wave))) = 0;
%   % show the filtered signal and interpolated signal
%     plot(t_canon',data_filtered, '.',t',data_interp,'-');

    %% Calculate PWM duty cycle
    period_samples = round(period * sample_rate);  % round舍入至最近的小数或整数                                  
    period_num = floor(N / period_samples); % floor向负无穷舍入
    duty = zeros(period_num, 1);
    pwm_wave = zeros(N, 1);
    for ii = 1:period_num
        strength = target_wave((ii-1) * period_samples + 1 );
        duty(ii) = (strength + 1) / 2 * (duty_upper_bound - duty_lower_bound) + duty_lower_bound;  % 将声音信号从-1~1映射到0~1，再映射到0.99~0.01之间
        busy_samples = round(duty(ii) * period_samples);
        pwm_wave((ii-1) * period_samples + 1: (ii-1) * period_samples + busy_samples) = 1;
    end
    duty = (duty-0.5)/max(abs(duty-0.5))*0.5+0.5;  %将duty扩大到0~1之间
    
    busy_time = duty * full_busy; % us
    idle_time = (1 - duty) * full_idle;
    results = [round(busy_time)];
    figure(4);
    plot (results);
    title(['results']);

    figure(5);
    plot (duty);
    title(['duty']);
    
%% 看信号频谱分布
figure(6);
n = length(results);
y0 = abs(fft(results,n)/n);
f3 = (0:n-1)*(target_frequency/n);
% y3 = abs(y0).^2/n;

% N2 = period_num;
% y0 = abs(fft(results)/N2);
% y3 = y0(1:N2/2+1);
% f = (0:round(N2/2))*target_frequency/N2;
plot(f3(2:length(results)/2),y0(2:length(results)/2))
% xlim=([-1, 3000])
title(['spectrogram of pwm wave']);



% N1 = time_len * fs;
% y1 = abs(fft(data)/N1);
%     y2 = y1(1:N1/2+1)
%     f1 = (0:round(N1/2))*fs/N1;
%     plot(f1,y2)
%     title('Spectrogram of AM modulated signal')



    %% Write duty cycle traces to txt files
    fid = fopen(['D:\OneDrive - zju.edu.cn\Research\Capspeaker期刊\聚川师兄资料\程序\调制程序-这几个是调制的程序，输入是低通滤波后的音频，输出是占空比\traces\',filename(1:end),'AM_20k_duty_cycle_32k.txt'], 'w');
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
    
end
toc