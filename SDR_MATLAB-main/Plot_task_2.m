% 清理环境
clear; clc; close all;

% 加载 Tracking 结果
load('TckResult_EphOpensky_90.mat', 'TckResult_Eph'); % 替换文件名

% 获取有效卫星索引（去除无数据的索引）
validSatellites = [3, 4, 16, 22, 26, 27, 31, 32]; % 只取这些卫星的数据
numSatellites = length(validSatellites);

%% **1. 绘制 Prompt 相关值 P_i**
figure;
hold on;
for i = 1:numSatellites
    satIndex = validSatellites(i); % 获取真实索引
    P_i_data = TckResult_Eph(satIndex).P_i;

    % 计算该卫星的时间轴，确保匹配数据长度
    timeAxis = 1:length(P_i_data);

    % 绘制该卫星数据
    plot(timeAxis, P_i_data, 'LineWidth', 1.5, 'DisplayName', ['Satellite ', num2str(satIndex)]);
end
xlabel('Time (ms)');
ylabel('Prompt Correlation (P_i)');
title('Prompt Correlation for Selected Satellites');
legend;
grid on;

%% **2. 绘制 P_q 相关值**
figure;
hold on;
for i = 1:numSatellites
    satIndex = validSatellites(i);
    P_q_data = TckResult_Eph(satIndex).P_q;
    timeAxis = 1:length(P_q_data);
    plot(timeAxis, P_q_data, 'LineWidth', 1.5, 'DisplayName', ['Satellite ', num2str(satIndex)]);
end
xlabel('Time (ms)');
ylabel('P_q Correlation');
title('P_q Correlation for Selected Satellites');
legend;
grid on;

%% **3. 绘制载波频率变化**
figure;
hold on;
for i = 1:numSatellites
    satIndex = validSatellites(i);
    carrierFreq_data = TckResult_Eph(satIndex).carrierFreq;
    timeAxis = 1:length(carrierFreq_data);
    plot(timeAxis, carrierFreq_data, 'LineWidth', 1.5, 'DisplayName', ['Satellite ', num2str(satIndex)]);
end
xlabel('Time (ms)');
ylabel('Carrier Frequency (Hz)');
title('Carrier Frequency Variation for Selected Satellites');
legend;
grid on;

disp('All plots have been generated successfully!');