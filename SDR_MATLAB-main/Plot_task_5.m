clc; clear; close all;

%% 1. 自动检测 Task 5 结果文件
baseFileName = 'Opensky'; % 替换为实际的 file.fileName
fileList = { ...
    ['navSolVT_', baseFileName, '_updated.mat'], ...
    ['navSolVT_mltCorr_', baseFileName, '_updated.mat'] ...
};

dataFileName = ''; % 初始化为空
for i = 1:length(fileList)
    if exist(fileList{i}, 'file')
        dataFileName = fileList{i};
        break;
    end
end

if isempty(dataFileName)
    error('未找到 Task 5 结果文件，请检查文件名或路径！');
end

disp(['加载文件: ', dataFileName]);
load(dataFileName, 'navSolutionsVT'); % 载入 Kalman 估计结果

%% 2. 提取关键数据
time = navSolutionsVT.localTime; % 时间
usrPosENU = navSolutionsVT.usrPosENU; % 用户位置 (ENU)
usrVel = navSolutionsVT.usrVel; % 用户速度 (ECEF)
usrVelENU = navSolutionsVT.usrVelENU; % 用户速度 (ENU)
clkDrift = navSolutionsVT.clkDrift; % 钟漂
clkBias = navSolutionsVT.clkBias; % 钟偏

%% 🎯 3. 绘制用户轨迹 (ENU 坐标系)
figure;
plot(usrPosENU(:,1), usrPosENU(:,2), 'b.-', 'LineWidth', 1.5);
xlabel('East (m)');
ylabel('North (m)');
title('User Trajectory in ENU Coordinates');
grid on;
axis equal;
legend('Trajectory');
saveas(gcf, 'Task5_User_Trajectory.png');

%% 🎯 4. 绘制速度变化 (XYZ 方向)
figure;
subplot(3,1,1);
plot(time, usrVel(:,1), 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Vx (m/s)');
title('Velocity in X Direction');

subplot(3,1,2);
plot(time, usrVel(:,2), 'g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Vy (m/s)');
title('Velocity in Y Direction');

subplot(3,1,3);
plot(time, usrVel(:,3), 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Vz (m/s)');
title('Velocity in Z Direction');

saveas(gcf, 'Task5_User_Velocity.png');

%% 🎯 5. 计算位置误差 (假设真实位置为均值)
truePos = mean(usrPosENU, 1); % 计算均值作为参考
posError = sqrt(sum((usrPosENU - truePos).^2, 2)); % 计算误差

figure;
plot(time, posError, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position Error (m)');
title('Position Error Over Time');
grid on;
saveas(gcf, 'Task5_Position_Error.png');

%% 🎯 6. 误差直方图（多路径效应分析）
figure;
histogram(posError, 50, 'FaceColor', 'b', 'FaceAlpha', 0.5);
xlabel('Position Error (m)');
ylabel('Frequency');
title('Position Error Distribution');
grid on;
saveas(gcf, 'Task5_Error_Histogram.png');

%% 🎯 7. 绘制钟漂 (Clock Drift) 变化
figure;
plot(time, clkDrift, 'g.-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Clock Drift (m/s)');
title('Clock Drift Over Time');
grid on;
legend('Clock Drift');
saveas(gcf, 'Task5_Clock_Drift.png');

%% 🎯 8. 绘制钟偏 (Clock Bias) 变化
figure;
plot(time, clkBias, 'm.-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Clock Bias (m)');
title('Clock Bias Over Time');
grid on;
legend('Clock Bias');
saveas(gcf, 'Task5_Clock_Bias.png');

disp('✅ Task 5 可视化完成，所有图像已生成！');