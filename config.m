function cfg = config()
% 3D 双雷达自校准仿真 - 信号级回波版本(唯一参数入口)
%
% 参数来源: Echo_Simulation_Design.md(式 E1-E13)
% 结构: cfg.room / cfg.radar(N站数组) / cfg.waveform / cfg.target / cfg.sim

    c = 299792458;                     % 光速 (m/s)

    %% 1. 房间参数
    cfg.room.length = 8;      % x 向长度 (m)
    cfg.room.width  = 6;      % y 向宽度 (m)
    cfg.room.height = 2.8;    % 天花板高 (m)

    %% 2. 波形参数 (FMCW, 式 E1)
    cfg.waveform.fc       = 60e9;     % 载频 (Hz)
    cfg.waveform.B        = 1.5e9;      % 带宽 (Hz)
    cfg.waveform.Tc       = 50e-6;    % chirp 时宽 (s)
    cfg.waveform.Fs       = 5e6;      % 快时间采样率 (Hz)
    cfg.waveform.N_chirps = 8;       % 每 TX chirp 数
    cfg.waveform.frame_rate = 50;      % 目标帧率 (Hz); 0 = 背靠背(帧周期 = 采集时长)
                                      % >0 时帧间自动插空闲: T_idle = max(0, 1/f - 采集时长)
    lambda = c / cfg.waveform.fc;

    %% 3. 雷达位姿与阵列参数 (N 站数组, 本期 N=2)
    % 真实位置: 天花板角内缩 0.3 m; 安装偏差: 站 1 全零, 站 2 有偏差
    p_true = [0.3, 0.3, 2.8;
              7.7, 0.3, 2.8];
    d_yaw   = [0, 3]';    % deg
    d_pitch = [0, 2]';    % deg
    d_roll  = [0, 0]';    % deg (预留)

    cfg.array.type = 'L'; % 'L': TX竖臂(测俯仰)+RX横臂(测方位), ref 基线
                          % 'ULA': TX/RX 共线横臂(测 u_y), 新配准算法入口用

    ntx = 2;              % TX 个数 (L 型竖臂, 沿站系 z 轴, 测俯仰)
    nrx = 4;              % RX 个数 (L 型横臂, 沿站系 y 轴, 测方位)
    N = size(p_true, 1);
    for i = 1:N
        cfg.radar(i).p_true  = p_true(i, :);
        cfg.radar(i).d_yaw   = d_yaw(i);
        cfg.radar(i).d_pitch = d_pitch(i);
        cfg.radar(i).d_roll  = d_roll(i);

        % L 型阵列: 位姿原点 = 阵列相位中心(阵中心), 式 E10
        cfg.radar(i).ntx = ntx;
        cfg.radar(i).nrx = nrx;
        cfg.radar(i).d_tx = lambda/2;   % TX 阵元间距 (m, 半波长)
        cfg.radar(i).d_rx = lambda/2;   % RX 阵元间距 (m, 半波长)

        % 单元方向图(式 E12): 锥角支撑 + cos^g_p 滚降, 标定到 -10 dB @ 60 deg
        cfg.radar(i).fov = 90;          % 方向图支撑半锥角 (deg)
        cfg.radar(i).g_p = 3.32;        % 滚降指数

        % 雷达方程常数(式 E7)
        cfg.radar(i).Pt_dbm = 12.5;     % 发射功率 (dBm)
        cfg.radar(i).G0_dbi = 8;        % 单元峰值增益 (dBi)
        cfg.radar(i).NF_db  = 15;       % 噪声系数 (dB)

        % 最大探测距离: 仅作波形校验与诊断, 非硬检测门限
        cfg.radar(i).max_range = 10;    % m
    end

    %% 4. 目标参数 (方体, 可见面散射中心, 式 E6)
    cfg.target.facet_size = 0.05;     % 面元尺寸 (m)
    cfg.target.phase_mode = 'mirror'; % 'mirror' 全零初相 / 'random' U[0,2pi) 种子生成
    cfg.target.targets(1).pos  = [3, 5, 0.1];   % 方体中心 (房间系, m)
    cfg.target.targets(1).size = [0.1, 0.1, 0.1];     % 三向尺寸 (m), 轴对齐
    % 朝向字段预留: cfg.target.targets(k).R(3,3), 本期固定 eye(3)

    %% 4b. 目标轨迹 (main_motion 轨迹遍历用, 静态 main 忽略)
    % 轨迹 = 静态点序列: 每帧目标冻结在轨迹点上, 按静态模型生成回波
    % (停-走冻结近似, 多普勒可忽略论证见 Echo_Simulation_Design.md 运动目标节)
    % 轨迹高度不单独设置, 自动取 targets(1).pos(3)(Scene_trajectory 联动)
    cfg.target.traj.type        = 'Z';    % 轨迹形状: 'Z' | 'C' | 'M'
    cfg.target.traj.speed       = 2.0;    % 匀速 (m/s)
    cfg.target.traj.n_frames_cap = 600;   % 帧数上限 (耗时控制)

    %% 5. 信号处理参数 (角度 FFT 补零)
    cfg.proc.nfft_az = 256;   % 方位维(rx 横臂)角度 FFT 补零长度
    cfg.proc.nfft_el = 256;   % 俯仰维(tx 竖臂)角度 FFT 补零长度

    %% 6. 仿真参数
    cfg.sim.seed = 42;         % 随机种子(可复现)
    cfg.sim.save = false;       % 是否存 cube 与 .mat 到 results/
    % cfg.sim.plot = false;      % 画图总开关: false 只算不画(轨迹遍历循环内恒不画图)
    cfg.sim.plot = true; 

    %% 7. 运动检测判定 (main_motion 用)
    cfg.motion.snr_th = 12;   % 峰 SNR (dB, 相对谱中值) 高于此判检出
end
