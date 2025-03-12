clear; 
clc;
format long g;

% 添加路径（如果需要）
addpath geo
addpath acqtckpos 

%% 参数初始化
[file, signal, acq, ~, ~, ~] = initParameters();

%% 执行 Acquisition
fprintf('Starting Acquisition...\n');
Acquired = acquisition(file, signal, acq);

%% 检查是否获取到了卫星
if isempty(Acquired.sv)
    fprintf('No satellites acquired. Check parameter settings.\n');
    return;
end

%% 绘制 Acquisition 结果
fprintf('Plotting Acquisition Results...\n');

% 创建图形窗口
figure;

% 子图1: SNR vs. SV ID
subplot(3,1,1);
bar(Acquired.sv, Acquired.SNR, 'b'); % SNR 条形图
xlabel('Satellite PRN');
ylabel('SNR (dB)');
title('Signal-to-Noise Ratio (SNR) of Acquired Satellites');
grid on;

% 子图2: Doppler Shift vs. SV ID
subplot(3,1,2);
stem(Acquired.sv, Acquired.Doppler, 'r', 'LineWidth', 1.5); % 多普勒频移
xlabel('Satellite PRN');
ylabel('Doppler Shift (Hz)');
title('Doppler Shift of Acquired Satellites');
grid on;

% 子图3: Code Delay vs. SV ID
subplot(3,1,3);
plot(Acquired.sv, Acquired.codedelay, 'go-', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Satellite PRN');
ylabel('Code Delay (samples)');
title('Code Delay of Acquired Satellites');
grid on;

% 总标题
sgtitle('GNSS Acquisition Results_Opensky');

fprintf('Acquisition and plotting completed.\n');