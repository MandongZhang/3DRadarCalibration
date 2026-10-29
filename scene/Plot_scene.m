function fig = Plot_scene(cfg, scene, tgt)
% 3D 场景图(tiledlayout 1x2): 左 = 房间/雷达/方体面元, 右 = L 型阵列放大
% 返回 figure 句柄供 main 导出存盘
%
% 面元按站可见性着色: 红 = 仅站 1 照明, 蓝 = 仅站 2, 品红 = 双站, 灰 = 无站。
% 可见性用照明近似: n·(p_radar - p_m) > 0(TX/RX 同位, 式 E13 的单站简化)。

    L = cfg.room.length;  W = cfg.room.width;  H = cfg.room.height;

    % 面元聚合与各站照明
    P   = cat(1, tgt.facet.pos).';      % (3,M) 房间系
    Nrm = cat(1, tgt.facet.nrm).';      % (3,M)
    nfac = size(P, 2);
    lit = false(nfac, numel(scene.radar));
    for i = 1:numel(scene.radar)
        lit(:, i) = sum(Nrm .* (scene.radar(i).p_true.' - P), 1).' > 0;
    end

    fig = figure('Name', '3D 场景布局', 'Color', 'w', 'Position', [50 50 1150 480]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, '天花板对角双雷达场景与方体目标');

    % ---------- 左: 3D 场景 ----------
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on'); box(ax1, 'on');

    % 房间: 地面 + 竖直棱 + 顶面
    plot3(ax1, [0 L L 0 0], [0 0 W W 0], [0 0 0 0 0], 'k-', 'LineWidth', 1);
    for cx = [0 L]
        for cy = [0 W]
            plot3(ax1, [cx cx], [cy cy], [0 H], 'k--');
        end
    end
    plot3(ax1, [0 L L 0 0], [0 0 W W 0], [H H H H H], 'k:', 'LineWidth', 0.5);

    % 方体线框
    for k = 1:numel(tgt)
        draw_box(ax1, tgt(k).pos - tgt(k).size/2, tgt(k).pos + tgt(k).size/2);
    end

    % 面元: 按照明组合着色
    col = zeros(nfac, 3);
    m1 = lit(:,1) & ~lit(:,2);   col(m1,:) = repmat([0.85 0.10 0.10], nnz(m1), 1);  % 红: 仅站 1
    m2 = ~lit(:,1) & lit(:,2);   col(m2,:) = repmat([0.10 0.30 0.85], nnz(m2), 1);  % 蓝: 仅站 2
    mb = lit(:,1) & lit(:,2);    col(mb,:) = repmat([0.70 0.10 0.70], nnz(mb), 1);  % 品红: 双站
    mn = ~lit(:,1) & ~lit(:,2);  col(mn,:) = repmat([0.6 0.6 0.6], nnz(mn), 1);     % 灰: 无站
    scatter3(ax1, P(1,:), P(2,:), P(3,:), 4, col, 'filled', ...
             'DisplayName', '面元(按照明着色)');

    % 雷达: 位置 + boresight + 视场锥边缘环
    cols = {'r', 'b'};
    for i = 1:numel(scene.radar)
        p = scene.radar(i).p_true;
        b = scene.radar(i).boresight;
        scatter3(ax1, p(1), p(2), p(3), 90, cols{i}, 'p', 'filled', ...
                 'DisplayName', sprintf('雷达 %d', i));
        quiver3(ax1, p(1), p(2), p(3), b(1)*2, b(2)*2, b(3)*2, 0, ...
                'Color', cols{i}, 'LineWidth', 1.6, 'HandleVisibility', 'off');
        th = linspace(0, 2*pi, 25);
        for a = [-scene.radar(i).fov, scene.radar(i).fov] * pi/180
            ring = [cos(th)*sind(a); sin(th)*sind(a); cos(a)*ones(size(th))] * 3;
            ring = scene.radar(i).R_true * ring + p.';
            plot3(ax1, ring(1,:), ring(2,:), ring(3,:), ':', 'Color', cols{i}, ...
                  'HandleVisibility', 'off');
        end
    end

    xlabel(ax1, 'x (m)'); ylabel(ax1, 'y (m)'); zlabel(ax1, 'z (m)');
    axis(ax1, 'equal'); xlim(ax1, [-0.5 L+0.5]); ylim(ax1, [-0.5 W+0.5]); zlim(ax1, [0 H+0.5]);
    view(ax1, 38, 22);
    legend(ax1, 'Location', 'northeastoutside', 'FontSize', 8);

    % ---------- 右: L 型阵列放大(站 1) ----------
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on'); axis(ax2, 'equal');
    rad = scene.radar(1);
    mm = 1e3;
    t_p = ((1:rad.ntx) - 1) * rad.d_tx * mm;                   % 站系 z (mm), 竖臂, 拐角起
    r_q = (1:rad.nrx) * rad.d_rx * mm;                         % 站系 y (mm), 横臂, 偏移一格
    plot(ax2, r_q, zeros(size(r_q)), 'bs', 'MarkerFaceColor', 'b', ...
         'MarkerSize', 7, 'DisplayName', sprintf('RX ×%d (沿站系 y, 横臂)', rad.nrx));
    plot(ax2, zeros(size(t_p)), t_p, 'ro', 'MarkerFaceColor', 'r', ...
         'MarkerSize', 7, 'DisplayName', sprintf('TX ×%d (沿站系 z, 竖臂)', rad.ntx));
    [VG, RG] = meshgrid(r_q, t_p);
    plot(ax2, VG(:), RG(:), 'x', 'Color', [0.5 0.5 0.5], ...
         'DisplayName', sprintf('TDM 虚拟阵 %dx%d', rad.ntx, rad.nrx));
    xlabel(ax2, '站系 y (mm)'); ylabel(ax2, '站系 z (mm)');
    title(ax2, sprintf('雷达 1 L 型阵列(阵元间距 %.2f mm)', rad.d_tx*mm));
    legend(ax2, 'Location', 'southeastoutside');
end

% ================================================================
function draw_box(ax, lo, hi)
% 轴对齐长方体 12 条棱
    x = lo(1);  X = hi(1);  y = lo(2);  Y = hi(2);  z = lo(3);  Z = hi(3);
    P = [x y z; X y z; X Y z; x Y z; x y Z; X y Z; X Y Z; x Y Z];
    seg = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 1 5; 2 6; 3 7; 4 8];
    for s = 1:size(seg, 1)
        plot3(ax, P(seg(s,:), 1), P(seg(s,:), 2), P(seg(s,:), 3), ...
              'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end
