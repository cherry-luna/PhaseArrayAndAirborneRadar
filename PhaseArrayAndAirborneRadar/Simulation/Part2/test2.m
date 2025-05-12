%% part2
clc;clear;close all;

%% 参数设定

c = 3e8;
eps = 0.0001;

% 雷达参数
N_R = 16;               % 阵元数
P_t = 2e3;              % 单阵元功率(W)
H_c = 5e3;              % 载机高度(m)
v_c = 150;              % 载机速度(m/s)
fc = 1e9;               % 载频(Hz)
lambda = c/fc;          % 波长(m)
B = 1e6;                % 带宽(Hz)
T_p = 100e-6;           % 脉宽(s)
PRF = 1e3;              % 脉冲重复频率(Hz)
CPI = 256;              % 积累脉冲数(CPI内脉冲数)
d = lambda/2;           % 阵元间距(m)
Ls_dB = 10;             % 接收机损耗(dB)
Ls = 10^(Ls_dB / 10);   % 转化为线性值
F_dB = 5;               % 噪声系数(dB)
F = 10^(F_dB / 10);     % 转化为线性值


% 目标参数
R_t = 90e3;        % 目标距离(m)
RCS_t = 5;         % 目标RCS(m²)
v_t = 60;          % 目标径向速度(m/s)

% 杂波参数
sigma0 = 0.01;          % 杂波后向散射系数
N_bin = 101;            % 杂波块个数
T0 = 290;               % 标准温度(K)
k_B = 1.38e-23;         % 玻尔兹曼常数

% 仿真参数
k_sam = 20;            % 样本个数（杂波距离环个数）
azimuth_target = 0;     % 目标方位（°）
azimuth_array = 90;     % 阵列与载机飞行方向夹角（°）

%% T5

% 获取杂波数据
[x_clutter, ~, ~, ~] = ClutterGen(H_c, R_t, v_c, azimuth_target, azimuth_array, ...
    N_bin, CPI, N_R, d, lambda, PRF, B, k_sam, sigma0, P_t, Ls);

% % 空时二维FFT
N_fft_s = 512; N_fft_d = 512;
% % 选择第一个距离环的数据
clutter_data = reshape(x_clutter(:,1), N_R, CPI);
st_spectrum = fftshift(fft2(clutter_data, N_fft_s, N_fft_d));
st_spectrum_db = 20 * log10(abs(st_spectrum.') + eps);

% 生成频率轴
fs_axis = linspace(-0.5, 0.5, N_fft_s);       % 归一化空间频率 (d/λ)
fd_axis = linspace(-0.5, 0.5, N_fft_d);       % 归一化多普勒频率 (Hz)

% 绘制三维空时谱
figure('Name','正侧视单距离环杂波空时谱', 'Position', [100, 100, 800, 500]);
surf(fs_axis, fd_axis, st_spectrum_db, 'EdgeColor', 'none', 'DisplayName','杂波脊');
xlabel('归一化空间频率 (d/\lambda)'); ylabel('归一化多普勒频率F_d (Hz)'); zlabel('幅度 (dB)');
% xlim([-0.5,0.5]);ylim([-0.5,0.5]);
title('正侧视单距离环杂波空时谱'); shading interp;colormap jet; colorbar;
axis tight;grid on;

% 理论杂波脊线叠加
hold on;
slope = (2 * v_c) / (d * PRF);% 理论斜率
fs_theo = linspace(-0.5, 0.5, 100);
fd_theo = slope * fs_theo ;
max_db = max(st_spectrum_db(:)); % 获取频谱最大幅度
plot3(fs_theo, fd_theo, max_db * ones(size(fs_theo)), 'r--', 'LineWidth', 2, 'DisplayName','理论脊');
xlim([-0.5,0.5]);ylim([-0.5,0.5]);
% legend('杂波谱', '理论杂波脊');
legend show;
hold off;

%% T6

% 信号波形生成
fs = 2 * B;
Ts = 1 / fs;
td = 2 * R_t / c; % 目标时延 (s)
N_d = round(td * fs); % 目标时延采样点数
fd = 2 * v_t / lambda; % 目标多普勒频率
K = B/T_p;                          % 调频斜率
t_chirp = linspace(0, T_p-Ts, T_p*fs); 
St = exp(1j*pi*K*t_chirp.^2);       % LFM信号
N_st = length(St);                  % 单个脉冲采样点数
N_PRT = round(1/PRF * fs);          % 单个PRT采样点数
num_pulses = 256;                   % 脉冲数

% 阵列接收信号生成（含目标、杂波、噪声）
array_phase = exp(-1j * 2 * pi * d / lambda * sind(azimuth_target) * (0:N_R-1)).'; % 导向矢量

% 生成目标回波 (调用阵列接收信号函数)
rx_target = rx_array_airborneradar(St, array_phase, N_R, num_pulses, ...
    fd, PRF, td, N_PRT, N_d, N_st, P_t, RCS_t, R_t, lambda, Ls);

% 生成杂波
clutter = ClutterGen(H_c, R_t, v_c, azimuth_target, azimuth_array, ...
    N_bin, num_pulses*N_PRT, N_R, d, lambda, PRF, B, k_sam, sigma0, P_t, Ls);
rx_clutter = reshape(clutter(:,1), num_pulses*N_PRT, N_R);

% 生成噪声
noise_power = F * k_B * T0 * B;
noise = sqrt(noise_power/2)*(randn(size(rx_target)) + 1j*randn(size(rx_target)));

% 合成接收信号
rx_array_signal = rx_target + rx_clutter + noise; % [N_PRT*num_pulses, N_R]

% DBF处理
weights_hamming = hamming(N_R);      % 汉明窗
steering_vector = array_phase .* weights_hamming; 
w = steering_vector / (array_phase' * steering_vector); % 归一化权值
rx_beamformed = rx_array_signal * conj(w); % 波束形成

% 分帧处理
rx_pulses = reshape(rx_beamformed, N_PRT, num_pulses);

% 脉冲压缩
% 雷达发射LFM信号后接收信号的功率（雷达方程计算）
Pr = (P_t * lambda^2 * RCS_t) / ((4*pi)^3 * R_t^4 * Ls);
Ar = sqrt(Pr);
h_mf = conj(fliplr(Ar * St)).'; % 匹配滤波器
dbf_mf_output = zeros(N_PRT + N_st - 1, num_pulses);
for i = 1:num_pulses
    dbf_mf_output(:,i) = conv(rx_pulses(:,i), h_mf, "full"); 
end

% 距离轴校正（补偿匹配滤波引入的延迟）
range_axis = ((0:size(dbf_mf_output,1)-1) - (N_st-1)) * c/(2*fs);
valid_idx = range_axis >= 0;
range_axis = range_axis(valid_idx) / 1e3; % 转换为km
dbf_mf_output = dbf_mf_output(valid_idx,:);

% MTD处理（多普勒FFT）
mtd_output = fftshift(fft(dbf_mf_output, [], 2), 2);

mtd_output_abs = abs(mtd_output);
mtd_output_db = 20*log10(mtd_output_abs / max(mtd_output_abs(:)) + eps);

% 速度轴计算
doppler_axis = (-num_pulses/2:num_pulses/2-1) * PRF / num_pulses; 
speed_axis = doppler_axis * lambda/2;

% 三维可视化
figure('Name','机载雷达检测目标接收信号三维图', 'Position', [100, 100, 800, 500]);
[X,Y] = meshgrid(range_axis, speed_axis);
surf(X, Y, mtd_output_db.', 'EdgeColor','none');
shading interp; colormap jet; colorbar;
xlabel('距离 (km)'); ylabel('速度 (m/s)'); zlabel('幅度 (dB)');
title('距离-速度-幅度三维图');

% 标记目标
hold on;
[~, idx] = max(mtd_output_db(:));
[row, col] = ind2sub(size(mtd_output_db), idx);
detected_range = range_axis(row);
detected_speed = speed_axis(col);
plot3(detected_range, detected_speed, mtd_output_db(row,col), ...
    'rp', 'MarkerSize', 5, 'LineWidth', 2, 'MarkerFaceColor', 'r');
text_str = sprintf('(%.1f km, %.1f m/s)', detected_range, detected_speed);
text(detected_range, detected_speed, mtd_output_db(row,col)+3,...
    text_str, 'FontSize', 12,  'Color','r', 'HorizontalAlignment', 'center');
legend('最终波形','检测目标');

%% CFAR检测（二维CA-CFAR）

% CFAR参数设置
guard = [5, 5];     % 距离/速度保护单元
train = [10, 10];   % 距离/速度参考单元
P_fa = 1e-6;        % 虚警概率

% 执行CFAR检测
[detection_map, target_info] = cfar_2d(mtd_output_db, guard, train, P_fa,...
                                      range_axis, speed_axis, R_t, v_t);

% 打印检测结果
fprintf('===== 目标检测结果 =====\n');
fprintf('真实目标位置: %.2f km, %.2f m/s\n', target_info.TrueRange, target_info.TrueSpeed);
fprintf('检测目标位置: %.2f km, %.2f m/s\n', target_info.DetectedRange, target_info.DetectedSpeed);
fprintf('距离相对误差: %.2f%%\n', target_info.RangeError);
fprintf('速度相对误差: %.2f%%\n', target_info.SpeedError);

% 绘制检测结果
figure('Position',[100 100 800 500]);
imagesc(range_axis, speed_axis, detection_map.');
xlabel('距离 (km)'); ylabel('速度 (m/s)'); title('CFAR检测结果');
colormap jet; colorbar;
hold on;

% 标记真实目标位置
plot(target_info.TrueRange, target_info.TrueSpeed, 'go',...
    'MarkerSize',10,'LineWidth',2,'DisplayName','真实目标');
text_str = sprintf('(%.2f km, %.2f m/s)', target_info.TrueRange, target_info.TrueSpeed);
text(target_info.TrueRange, target_info.TrueSpeed+3, ...
    text_str, 'FontSize', 12, 'Color', 'g', 'HorizontalAlignment', 'center');

% 标记检测目标位置
if ~isnan(target_info.DetectedRange)
    plot(target_info.DetectedRange, target_info.DetectedSpeed, 'rp',...
        'MarkerSize',12,'LineWidth',2,'DisplayName','检测目标');
    text_str = sprintf('(%.2f km, %.2f m/s)', target_info.DetectedRange, target_info.DetectedSpeed);
    text(target_info.DetectedRange, target_info.DetectedSpeed-3, ...
    text_str, 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
    legend show;
end

%% T7

% 参数设置
v_max = PRF * lambda / 4; % 归一化多普勒频率±0.5对应的最大速度
v_values = linspace(-v_max, v_max, 31); % 生成速度范围(-v_max到v_max,共31个点)
SNCR_dB = zeros(size(v_values)); % 存储各速度对应的SCNR

% 背景功率计算 (无目标时)
rx_background = rx_clutter + noise;
rx_beamformed_bg = rx_background * conj(w);
rx_pulses_bg = reshape(rx_beamformed_bg, N_PRT, num_pulses);

% 脉冲压缩和MTD处理
dbf_mf_output_bg = zeros(N_PRT + N_st - 1, num_pulses);
for i = 1:num_pulses
    dbf_mf_output_bg(:,i) = conv(rx_pulses_bg(:,i), h_mf, "full"); 
end

% 距离轴校正（补偿匹配滤波引入的延迟）
range_axis_bg = ((0:size(dbf_mf_output_bg,1)-1) - (N_st-1)) * c/(2*fs);
valid_idx_bg = range_axis_bg >= 0;
range_axis_bg = range_axis_bg(valid_idx_bg) / 1e3; % 转换为km
dbf_mf_output_bg = dbf_mf_output_bg(valid_idx_bg,:);

% MTD处理（多普勒FFT）
mtd_output_bg = fftshift(fft(dbf_mf_output_bg, [], 2), 2);
mtd_output_bg_abs = abs(mtd_output_bg);
mtd_power_bg = mean(mtd_output_bg_abs(:).^2); % 背景平均功率

% 遍历目标速度计算SCNR
for k = 1:length(v_values)
    v_k = v_values(k);
    fd_k = 2 * v_k / lambda; % 当前速度对应的多普勒频率
    
    % 生成目标信号
    rx_target_k = rx_array_airborneradar(St, array_phase, N_R, num_pulses, ...
        fd_k, PRF, td, N_PRT, N_d, N_st, P_t, RCS_t, R_t, lambda, Ls);
    
    % 合成接收信号
    rx_array_signal_k = rx_target_k + rx_background;
    
    % 处理链
    rx_beamformed_k = rx_array_signal_k * conj(w);
    rx_pulses_k = reshape(rx_beamformed_k, N_PRT, num_pulses);
    
    % 脉冲压缩
    dbf_mf_output_k = zeros(N_PRT + N_st -1, num_pulses);
    for i = 1:num_pulses
        dbf_mf_output_k(:,i) = conv(rx_pulses_k(:,i), h_mf, "full"); 
    end
    dbf_mf_output_k = dbf_mf_output_k(valid_idx_bg,:);
    
    % MTD处理
    mtd_output_k = fftshift(fft(dbf_mf_output_k, [], 2), 2);
    mtd_output_k_abs = abs(mtd_output_k);
    mtd_power_k = mtd_output_k_abs.^2;
    
    % 定位目标单元
    [~, range_idx] = min(abs(range_axis_bg * 1e3 - R_t));
    [~, doppler_idx] = min(abs(speed_axis - v_k));
    
    % 计算SCNR
    signal_power = mtd_power_k(range_idx, doppler_idx);
    SNCR_dB(k) = 10*log10(signal_power / mtd_power_bg);
end

% 绘制SCNR曲线
figure('Position',[100 100 800 500]);
plot(v_values, SNCR_dB, 'LineWidth',1.5);
xlabel('目标速度 (m/s)'); ylabel('SCNR (dB)');
title('目标SCNR随速度变化曲线');
grid on;
hold on;

% 标记杂波脊位置
v_clutter_ridge = v_c * sind(azimuth_target); % 正侧视杂波速度
plot([v_clutter_ridge, v_clutter_ridge], ylim, 'r--', 'LineWidth',1.2);
legend('SCNR曲线', '杂波脊位置');
