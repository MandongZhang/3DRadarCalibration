function g = Antenna_gain(u, boresight, fov_deg, g_p)
% 单元功率方向图(式 E12): g = cos^g_p(离轴角), 视场(锥角)外为 0
%
% 输入: u (3,M)        面元方向单位向量(非归一化时内部归一)
%       boresight (3,1) 阵列法向(站系 -> 房间系后的波束指向)
%       fov_deg         方向图支撑半锥角 (deg)
%       g_p             滚降指数(标定到 -10 dB @ 60 deg 时 g_p = 3.32)
% 输出: g (M,1) 线性功率增益(峰值 1), 视场外为 0
%
% 旧式(8)的 cos^g_p 单元滚降在信号级作用于逐面元逐阵元, 而非量测级近似

    u = u ./ max(vecnorm(u, 2, 1), eps);       % 保险归一 (3,M)
    cp = u.' * boresight(:);                   % (M,1)
    psi = acosd(min(max(cp, -1), 1));          % 离轴角 (deg)
    g = max(cosd(psi), 0).^g_p;
    g(psi > fov_deg) = 0;
end
