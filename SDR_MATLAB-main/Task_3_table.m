%% Ephemeris Table Export from .MAT Files
clear; clc;

% 选择.mat文件
[filename, pathname] = uigetfile('*.mat', 'Select Ephemeris Data File');
if isequal(filename, 0)
    error('No file selected');
end
matFile = fullfile(pathname, filename);

% 加载MAT文件
try
    load(matFile, 'TckResult_Eph'); % 强制加载名为eph的变量
catch
    error('Failed to load eph structure from selected file');
end

% 验证数据结构
if ~exist('eph', 'var') || ~isstruct(eph)
    error('Invalid ephemeris structure in selected file');
end

% 用户输入PRN编号
prn = input('Enter target PRN number (e.g. 16): ');
if ~ismember(prn, [eph.PRN])
    error('PRN %d not found in the ephemeris data', prn);
end

%% 创建参数列表 (全英文)
parameters = {
    'GPS Week',         'GPS week number',                                  eph(prn).weekNumber;
    'Toe',              'Time of Ephemeris (seconds)',                      eph(prn).toe;
    'IODC/IODE',        'Issue of Data (Clock/Ephemeris)',                  [eph(prn).IODC, eph(prn).IODE];
    'Semi-Major Axis',  'Square root of semi-major axis (√a, m^0.5)',       eph(prn).sqrtA;
    'Eccentricity',     'Orbital eccentricity',                             eph(prn).e;
    'Inclination',      'Inclination angle at reference time (rad)',        eph(prn).i0;
    'RAAN',             'Right Ascension of Ascending Node (rad)',          eph(prn).omega0;
    'Argument of Perigee','Argument of perigee (rad)',                      eph(prn).omega;
    'Mean Anomaly',     'Mean anomaly at reference time (rad)',             eph(prn).M0;
    'Delta_n',          'Mean motion difference (rad/s)',                   eph(prn).delta_n;
    'RAAN Rate',        'Rate of RAAN (rad/s)',                             eph(prn).OMEGAdot;
    'ArgP Rate',        'Rate of argument of perigee (rad/s)',              eph(prn).omegaDot;
    'Cuc',              'Latitude argument cosine correction (rad)',        eph(prn).Cuc;
    'Cus',              'Latitude argument sine correction (rad)',          eph(prn).Cus;
    'TGD',              'Time group delay (seconds)',                       eph(prn).TGD;
};

%% 数值格式化
formatted_values = cell(size(parameters,1), 1);
for i = 1:size(parameters,1)
    try
        current_value = parameters{i,3};
        
        if contains(parameters{i,1}, 'IODC/IODE')
            formatted_values{i} = sprintf('%d/%d', current_value(1), current_value(2));
        elseif abs(current_value) < 1e-3
            formatted_values{i} = sprintf('%.4e', current_value);
        else
            formatted_values{i} = sprintf('%.16g', current_value);
        end
    catch
        formatted_values{i} = 'N/A';
    end
end

%% 生成表格
ephemerisTable = table(...
    parameters(:,1), ...
    parameters(:,2), ...
    formatted_values, ...
    'VariableNames', {'Parameter', 'Description', 'Value'});

% 显示表格
fprintf('\n\nTable 1: Ephemeris Parameters for PRN %d\n', prn);
disp(ephemerisTable);

%% 导出文件
[~, name] = fileparts(filename);
exportName = sprintf('%s_PRN%02d.csv', name, prn);
writetable(ephemerisTable, exportName);

fprintf('\nExport successful: %s\n', fullfile(pwd, exportName));