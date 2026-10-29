function tgt = Scene_targets(cfg)
% 方体目标 6 面面元网格(第一块: 场景布局)
%
% 输入: cfg - config.m 输出
% 输出: tgt(k) - 目标结构体数组
%   tgt(k).pos  (1,3)   方体中心 (房间系, m)
%   tgt(k).size (1,3)   三向尺寸 (m)
%   tgt(k).R    (3,3)   朝向(本期轴对齐 eye(3), 字段预留)
%   tgt(k).facet.pos   (M_k,3)  面元中心 (房间系, m), 位于面平面上
%   tgt(k).facet.nrm   (M_k,3)  面元外法向 (单位)
%   tgt(k).facet.area  (M_k,1)  面元实际面积 (m^2)
%   tgt(k).facet.phase (M_k,1)  散射初相: mirror 全 0 / random ~ U[0,2pi)
%   tgt(k).M                    该目标面元数
%
% 网格规则(Echo_Simulation_Design.md §3.1): 每边 n = round(L/facet_size)
% 等分, d_eff = L/n 精确铺满不留半格。面元尺寸相对波长(约 3.9 mm)大得多,
% 每个 面元视作小平板散射体, RCS 由式 E6 在 Echo_generate 中计算。
%
% 散射相位是目标属性: 全局 rng(由调用方设置) 一次生成, 不逐站重新生成。
% 面元对雷达的可见性(TX 照明 + RX 可见, 式 E13)留给 Echo_generate 逐站判断。

    % 判断是否设置了目标的尺寸大小
    if ~isfield(cfg.target, 'facet_size'), error('cfg.target.facet_size 缺失'); end


    for k = 1:numel(cfg.target.targets)
        tk = cfg.target.targets(k);
        C  = tk.pos(:).';
        L  = tk.size(:).';

        % 没有定义旋转，则默认朝向与房间坐标正交
        if isfield(tk, 'R') && ~isempty(tk.R), Rk = tk.R; else, Rk = eye(3); end

        % 6 个面: 法向 = ±三轴(Rk 各列)
        dirs = [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];   % (3,6) 轴对齐法向
        P = zeros(0,3);  Nrm = zeros(0,3);  A = zeros(0,1);
        for f = 1:6
            n  = dirs(:, f);                     % 面法向(轴对齐单位向量)
            nz = find(n ~= 0, 1);
            
            % 面内两轴: 其余两个轴向
            tang = setdiff(1:3, nz);
            t1 = zeros(3,1);  t1(tang(1)) = 1;   % 面内轴 1(单位轴向量; 修复原 dirs(:,tang) 用轴序号错索引 6 列方向表的 bug)
            t2 = zeros(3,1);  t2(tang(2)) = 1;   % 面内轴 2

            L1 = L(tang(1));  L2 = L(tang(2));
            n1 = max(round(L1 / cfg.target.facet_size), 1);
            n2 = max(round(L2 / cfg.target.facet_size), 1);
            d1 = L1 / n1;  d2 = L2 / n2;         % 精确铺满的格距

            u1 = ((1:n1) - 0.5) * d1 - L1/2;     % 格中心 (局部)
            u2 = ((1:n2) - 0.5) * d2 - L2/2;
            [U1, U2] = meshgrid(u1, u2);

            % 局部面坐标 -> 房间系: 中心 + 法向半尺寸 + R*(t1*u1 + t2*u2)
            loc = (t1 * U1(:).' + t2 * U2(:).');             % (3, n1*n2)
            Pf  = C.' + Rk * (n * L(nz)/2) + Rk * loc;
            Nf  = repmat((Rk * n).', n1*n2, 1);              % (n1*n2, 3) 外法向
            Af  = d1 * d2 * ones(n1*n2, 1);

            P   = [P;   Pf.'];                               %#ok<AGROW>
            Nrm = [Nrm; Nf];                                 %#ok<AGROW>
            A   = [A;   Af];                                 %#ok<AGROW>
        end

        tgt(k).pos  = C;
        tgt(k).size = L;
        tgt(k).R    = Rk;
        tgt(k).facet.pos   = P;
        tgt(k).facet.nrm   = Nrm;
        tgt(k).facet.area  = A;
        switch cfg.target.phase_mode
            case 'mirror'
                tgt(k).facet.phase = zeros(size(A));
            case 'random'
                tgt(k).facet.phase = 2*pi * rand(size(A));
            otherwise
                error('未知 phase_mode: %s', cfg.target.phase_mode);
        end
        tgt(k).M = size(P, 1);
    end
end
