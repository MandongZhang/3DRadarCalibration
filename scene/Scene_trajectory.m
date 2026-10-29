function traj = Scene_trajectory(cfg, wf)
% 目标轨迹生成(Z/C/M 三种形状, 弧长匀速采样到帧网格)
%
% 输入: cfg - config.m 输出
%       wf  - Waveform_derive 输出(帧周期 wf.T_frame 已含 cfg.waveform.frame_rate
%             决定的帧间空闲)
% 输出: traj
%   traj.type  形状类型 ('Z'|'C'|'M')
%   traj.t     (Nf,1)  帧时刻 (s),  t_n = (n-1)*T_frame
%   traj.pos   (Nf,3)  方体中心房间系坐标 (m), 每帧一个静态点(轨迹遍历用)
%   traj.wp    (Wp,3)  形状锚点(绘图参考)
%
% 帧周期: wf.T_frame = ntx*N_chirps*Tc + T_idle(采集时长 + 帧间空闲, 由目标帧率
% cfg.waveform.frame_rate 控制)。轨迹遍历把轨迹看成多个静态点: 每帧目标冻结在
% 帧起始位置, 按静态模型生成回波(停-走冻结)。近似成立条件与多普勒可忽略论证
% 见 Echo_Simulation_Design.md 运动目标节。
%
% 锚点默认值(房间系 xy, 高度取 cfg.target.targets(1).pos(3) 联动, 均在 8x6 房间内):
%   Z: [2.5 1.5] -> [5.5 1.5] -> [2.5 4.5] -> [5.5 4.5]   (三段折线)
%   C: 圆心 [4 3], r = 1.5, 弧 60 deg -> 300 deg          (开口朝 +x)
%   M: [1.8 1.8] -> [2.9 4.2] -> [4.0 1.8] -> [5.1 4.2] -> [6.2 1.8]  (四段折线)

    tr = cfg.target.traj;
    T_frame = wf.T_frame;

    % ---- 形状锚点(房间系 xy) ----
    switch upper(tr.type)
        case 'Z'
            wp = [2.5 1.5; 5.5 1.5; 2.5 4.5; 5.5 4.5];
        case 'C'
            th = linspace(60, 300, 81).';                  % deg, 密集采样近似圆弧
            wp = [4, 3] + 1.5 * [cosd(th), sind(th)];
        case 'M'
            wp = [1.8 1.8; 2.9 4.2; 4.0 1.8; 5.1 4.2; 6.2 1.8];
        otherwise
            error('未知轨迹类型: %s (可选 Z/C/M)', char(tr.type));
    end
    z0 = cfg.target.targets(1).pos(3);                     % 轨迹高度 = 目标1定义高度(自动联动)
    Wp = [wp, z0 * ones(size(wp, 1), 1)];                  % (Wp,3) 补高度

    % ---- 折线弧长 ----
    seg   = vecnorm(diff(Wp, 1, 1), 2, 2);                 % (Wp-1,1) 各段长度
    s_cum = [0; cumsum(seg)];                              % (Wp,1) 累计弧长
    L_path = s_cum(end);

    % ---- 匀速采样到帧网格(轨迹末端之后冻结在终点) ----
    T_total = L_path / tr.speed;
    Nf = min(floor(T_total / T_frame) + 1, tr.n_frames_cap);
    traj.t  = (0:Nf-1).' * T_frame;                        % (Nf,1)
    s_n     = min(tr.speed * traj.t, L_path);              % 每帧的弧长坐标
    traj.pos(:, 1) = interp1(s_cum, Wp(:, 1), s_n, 'linear');
    traj.pos(:, 2) = interp1(s_cum, Wp(:, 2), s_n, 'linear');
    traj.pos(:, 3) = interp1(s_cum, Wp(:, 3), s_n, 'linear');

    traj.type = upper(tr.type);
    traj.wp   = Wp;
end
