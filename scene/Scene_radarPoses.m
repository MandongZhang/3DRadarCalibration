 function scene = Scene_radarPoses(cfg)
% 双雷达位姿 + L 型阵列几何(第一块: 场景布局)
%
% 输入: cfg - config.m 输出
% 输出: scene - 场景结构体
%   scene.radar(i).p_true / R_true     真实位姿(量测用真值)
%   scene.radar(i).p_prior / R_prior   理想对角先验位姿(偏差置零, 二期求解初值)
%   scene.radar(i).tx_pos (3,ntx)      TX 阵元房间系坐标 = p_true + R*[0; t_p; 0]
%   scene.radar(i).rx_pos (3,nrx)      RX 阵元房间系坐标 = p_true + R*[0; 0; r_q]
%   scene.radar(i).boresight (3,1)     = R(:,1), 波束指向
%   scene.radar(i).fov / g_p / Pt_dbm / G0_dbi / NF_db / max_range  透传
%   scene.radar(i).ntx / nrx / d_tx / d_rx                          透传
%
% 坐标系约定(Echo_Simulation_Design.md §0):
%   房间系 W: 原点墙角, x 长 y 宽 z 向上
%   站系 S:   前 x(波束指向) 右 y 上 z
%   C^s = R' * (C - p),  R 为站系到房间系旋转
%
% 旋转约定: R = Rz(yaw) * Ry(pitch) * Rx(roll), pitch > 0 为下视
%
% L 型阵列(式 E10): 位姿原点 = L 型拐角(TX1 位置)
%   TX 沿站系 z 轴正向(竖臂, 测俯仰): t_p = (p - 1) * d_tx
%   RX 沿站系 y 轴正向(横臂, 测方位): r_q = q * d_rx(整体偏移一个阵元, 拐角处
%   TX1 与 RX1 间隔 d_rx, 避免收发阵元位置重叠)
%   TDM-MIMO 虚拟阵元 = t_p + r_q(Echo_generate 元数据)

    for i = 1:numel(cfg.radar)
        p = cfg.radar(i).p_true;

        % 理想朝向: 指向地面中心 [L/2, W/2, 0]
        d = [cfg.room.length/2 - p(1), cfg.room.width/2 - p(2), -p(3)];
        ideal_yaw   = atan2d(d(2), d(1));          % 水平方位
        ideal_pitch = -asind(d(3) / norm(d));      % 下视角为正(见旋转约定)

        % 真实位姿 = 理想 + 安装偏差
        yaw_t   = ideal_yaw   + cfg.radar(i).d_yaw;
        pitch_t = ideal_pitch + cfg.radar(i).d_pitch;
        roll_t  = cfg.radar(i).d_roll;

        scene.radar(i).p_true = p;
        scene.radar(i).R_true = euler2R(yaw_t, pitch_t, roll_t);

        % 先验位姿: 偏差置零
        scene.radar(i).p_prior = p;
        scene.radar(i).R_prior = euler2R(ideal_yaw, ideal_pitch, 0);

        % 阵元坐标(站系 -> 房间系, 式 E10): 位姿原点 = 阵列参考元(TX1)
        ntx = cfg.radar(i).ntx;   nrx = cfg.radar(i).nrx;
        if strcmp(cfg.array.type, 'ULA')
            % ULA: TX/RX 共线沿站系 y; TX 长距间隔 nrx*d_rx, RX 密集 d_rx
            % 虚拟阵元 = t_p + r_q, 恰为 0..(ntx*nrx-1)*d 连续半波长栅格
            % (TX1 与 RX1 位置重合在 TDM 分时下允许)
            t_p = (0:ntx-1) * nrx * cfg.radar(i).d_rx;         % (1,ntx) 沿 y
            r_q = (0:nrx-1) * cfg.radar(i).d_rx;               % (1,nrx) 沿 y, 从 0 起
            tx_local = [zeros(1,ntx); t_p; zeros(1,ntx)];
            rx_local = [zeros(1,nrx); r_q; zeros(1,nrx)];
        else
            % L 型(ref): TX 竖臂沿 z 测俯仰, RX 横臂沿 y 测方位; RX 臂偏移一格
            t_p = ((1:ntx) - 1) * cfg.radar(i).d_tx;           % (1,ntx) 沿 z(竖臂), 从 0 起
            r_q = (1:nrx) * cfg.radar(i).d_rx;                 % (1,nrx) 沿 y(横臂), 从 d_rx 起
            tx_local = [zeros(1,ntx); zeros(1,ntx); t_p];      % (3,ntx) 站系
            rx_local = [zeros(1,nrx); r_q; zeros(1,nrx)];      % (3,nrx) 站系
        end
        scene.radar(i).tx_pos = p.' + scene.radar(i).R_true * tx_local;
        scene.radar(i).rx_pos = p.' + scene.radar(i).R_true * rx_local;
        scene.radar(i).boresight = scene.radar(i).R_true(:, 1);

        % 参数透传
        scene.radar(i).ntx = ntx;      scene.radar(i).nrx = nrx;
        scene.radar(i).d_tx = cfg.radar(i).d_tx;
        scene.radar(i).d_rx = cfg.radar(i).d_rx;
        scene.radar(i).fov  = cfg.radar(i).fov;
        scene.radar(i).g_p  = cfg.radar(i).g_p;
        scene.radar(i).Pt_dbm = cfg.radar(i).Pt_dbm;
        scene.radar(i).G0_dbi = cfg.radar(i).G0_dbi;
        scene.radar(i).NF_db  = cfg.radar(i).NF_db;
        scene.radar(i).max_range = cfg.radar(i).max_range;
    end
end

function R = euler2R(yaw, pitch, roll)
% ZYX 欧拉角 -> 旋转矩阵(角度制)
    cy = cosd(yaw);   sy = sind(yaw);
    cp = cosd(pitch); sp = sind(pitch);
    cr = cosd(roll);  sr = sind(roll);
    Rz = [cy -sy 0; sy cy 0; 0 0 1];
    Ry = [cp 0 sp; 0 1 0; -sp 0 cp];
    Rx = [1 0 0; 0 cr -sr; 0 sr cr];
    R = Rz * Ry * Rx;
end
