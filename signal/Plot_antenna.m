function fig = Plot_antenna(cfg)
% 单元方向增益图(诊断): config 的 fov / g_p / G0_dbi 可视化
%   左: 当前参数的 dB 滚降曲线, 标 fov 支撑边界与 -10 dB @ 60 deg 标定点
%   右: 不同滚降指数对比(线性), 标注各自 -3 dB 半功率角
% 返回 figure 句柄供 main 导出存盘

    rad = cfg.radar(1);
    psi = 0:0.1:90;                                  % (1,901) 行向量
    u = [cosd(psi); sind(psi); zeros(size(psi))];    % (3,901) 相对 boresight [1;0;0] 的离轴扫描

    fig = figure('Name', '单元方向图', 'Color', 'w', 'Position', [80 80 980 420]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('单元功率方向图 G = cos^{%.2f}\\psi(峰值 %d dBi)', ...
          rad.g_p, rad.G0_dbi));

    % ---- 左: 当前参数 dB 滚降 ----
    ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
    g1 = Antenna_gain(u, [1; 0; 0], rad.fov, rad.g_p);
    plot(ax1, psi, 10*log10(max(g1, 1e-12)), 'b-', 'LineWidth', 1.8, ...
         'DisplayName', sprintf('g_p = %.2f', rad.g_p));
    yline(ax1, -3, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(ax1, -10, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(ax1, rad.fov, 'r--', 'LineWidth', 1.3, ...
          'DisplayName', sprintf('fov = %d deg(支撑外置 0)', rad.fov));
    plot(ax1, 60, -10, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 7, ...
         'DisplayName', '-10 dB @ 60 deg 标定点');
    text(ax1, rad.fov + 2, -20, sprintf('fov=%d°', rad.fov), 'Color', 'r');
    ylim(ax1, [-40 2]); xlim(ax1, [0 90]);
    xlabel(ax1, '离轴角 \psi (deg)'); ylabel(ax1, 'G (dB, 归一)');
    title(ax1, 'dB 域滚降(config 当前参数)');
    legend(ax1, 'Location', 'southwest');

    % ---- 右: 滚降指数对比(线性) ----
    ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
    cols = lines(3);
    gps = [2, rad.g_p, 5];
    for k = 1:numel(gps)
        gk = Antenna_gain(u, [1; 0; 0], rad.fov, gps(k));
        plot(ax2, psi, gk, 'LineWidth', 1.5, 'Color', cols(k, :), ...
             'DisplayName', sprintf('g_p = %.2f(-3dB @ %.1f°)', ...
             gps(k), acosd(0.5^(1/gps(k)))));
    end
    plot(ax2, [0 90], [0.5 0.5], 'k:', 'HandleVisibility', 'off');
    ylim(ax2, [0 1.05]); xlim(ax2, [0 90]);
    xlabel(ax2, '离轴角 \psi (deg)'); ylabel(ax2, 'G(归一)');
    title(ax2, '滚降指数对比(线性)');
    legend(ax2, 'Location', 'northeast');
end
