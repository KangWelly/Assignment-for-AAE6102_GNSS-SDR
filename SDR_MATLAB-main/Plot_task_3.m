clear; clc; close all;

% Load data
load('eph_Opensky_90.mat'); % Ephemeris data
load('TckResult_EphOpensky_90.mat'); % Prompt signal data

ephemeris = eph; % Ensure the variable is accessible
num_sv = length(ephemeris); % Get the number of satellites

figure;
hold on;

% --- 1. Plot eccentricity vs. Toe ---
subplot(2,2,1);
hold on;
for prn = 1:num_sv
    if isempty(ephemeris(prn).toe) || isempty(ephemeris(prn).ecc), continue; end
    plot(ephemeris(prn).toe, ephemeris(prn).ecc, '-o', 'DisplayName', ['SV' num2str(prn)]);
end
xlabel('Toe (s)');
ylabel('Eccentricity (e)');
title('Satellite Orbital Eccentricity vs. Toe');
legend;
grid on;

% --- 2. Plot semi-major axis vs. Toe ---
subplot(2,2,2);
hold on;
for prn = 1:num_sv
    if isempty(ephemeris(prn).toe) || isempty(ephemeris(prn).sqrta), continue; end
    semi_major_axis = (ephemeris(prn).sqrta).^2;  % Semi-major axis: a = (sqrt(A))^2
    plot(ephemeris(prn).toe, semi_major_axis, '-o', 'DisplayName', ['SV' num2str(prn)]);
end
xlabel('Toe (s)');
ylabel('Semi-Major Axis (m)');
title('Semi-Major Axis vs. Toe');
legend;
grid on;

% --- 3. Plot inclination angle (i0) vs. Toe ---
subplot(2,2,3);
hold on;
for prn = 1:num_sv
    if isempty(ephemeris(prn).toe) || isempty(ephemeris(prn).i0), continue; end
    % Ensure matching dimensions for Toe and i0
    min_len = min(length(ephemeris(prn).toe), length(ephemeris(prn).i0));
    plot(ephemeris(prn).toe(1:min_len), rad2deg(ephemeris(prn).i0(1:min_len)), '-o', 'DisplayName', ['SV' num2str(prn)]);
end
xlabel('Toe (s)');
ylabel('Inclination Angle (°)');
title('Orbital Inclination Angle vs. Toe');
legend;
grid on;

% --- 4. Plot navigation data Prompt signal ---
subplot(2,2,4);
hold on;
for prn = 1:num_sv
    if prn <= length(TckResult_Eph) && isfield(TckResult_Eph, 'P_i')
        plot(TckResult_Eph(prn).P_i, 'DisplayName', ['SV' num2str(prn)]);
    end
end
xlabel('Samples');
ylabel('Prompt Signal Amplitude');
title('Navigation Data Prompt Signal');
legend;
grid on;

hold off;
sgtitle('Task 3: Navigation Data Decoding Results');