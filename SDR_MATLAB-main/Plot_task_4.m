%% 清理环境
clear; clc; close all;

%% 加载 Task 4 相关数据
fileName = 'Opensky'; % 这里的 Opensky 对应文件名
load(['navSolCT_10ms_', fileName, '.mat']); % 载入 WLS 位置解
load(['tckRstCT_10ms_', fileName, '.mat']); % 载入跟踪结果

%% 真实值（WGS-84 坐标系）
true_lat = 22.328444770087565;
true_lon = 114.1713630049711;
true_alt = 1; % 假设海拔高度为 0

% 转换真实值到 ECEF（地心地固坐标系）
true_xyz = lla2ecef([true_lat, true_lon, true_alt]); % 直接用 MATLAB 内置函数

%% 确保 navSolutionsCT 存在
if ~exist('navSolutionsCT', 'var') || isempty(navSolutionsCT)
error('Error: navSolutionsCT 数据为空，无法绘图。');
end

%% 提取位置和时间数据
time = (1:length(navSolutionsCT.usrPos)) * 10; % 10ms 采样时间步长
pos_ECEF = navSolutionsCT.usrPos; % 用户 ECEF 坐标 (X, Y, Z)
clkBias = navSolutionsCT.clkBias; % 钟差

%% 计算位置误差
pos_error = sqrt(sum((pos_ECEF - true_xyz).^2, 2)); % 计算位置误差 (m)

%% 计算速度（简单差分方法）
velocity = diff(pos_ECEF) ./ diff(time'); % 计算ECEF速度
velocity = [velocity; velocity(end, :)]; % 维持相同行数

%% 绘制用户 ECEF 轨迹
figure;
plot3(pos_ECEF(:,1), pos_ECEF(:,2), pos_ECEF(:,3), 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8); % 蓝色轨迹，点线结合
hold on;
scatter3(true_xyz(1), true_xyz(2), true_xyz(3), 100, 'r', 'filled'); % 真实值
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('User Trajectory (ECEF)');
legend('WLS Estimated Trajectory', 'Ground Truth', 'Location', 'best');
axis tight; % 自动调整坐标轴范围
hold off;

%% 绘制速度变化 (XYZ 方向)
figure;
subplot(3,1,1);
plot(time, velocity(:,1), 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Vx (m/s)');
title('Velocity in X Direction');

subplot(3,1,2);
plot(time, velocity(:,2), 'g', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Vy (m/s)');
title('Velocity in Y Direction');

subplot(3,1,3);
plot(time, velocity(:,3), 'r', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Vz (m/s)');
title('Velocity in Z Direction');

%% 位置误差分析
figure;
plot(time, pos_error, 'k', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Position Error (m)');
title('Position Error Over Time');
grid on;

%% 钟差变化
figure;
plot(time, clkBias, 'b', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Clock Bias (m)');
title('Clock Bias Over Time');
grid on;

%% 误差直方图（多路径效应分析）
figure;
histogram(pos_error, 50);
xlabel('Position Error (m)');
ylabel('Frequency');
title('Position Error Distribution (Multipath Effect Analysis)');
grid on;

disp('Task 4 Result Visualization Completed!');