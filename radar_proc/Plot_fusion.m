function fig = Plot_fusion(traj, P_hist, T_hist, both, fuse, R1, p1)
% 双站点云融合结果图(tiledlayout 1x2, 构造房间系视角)
%   左: 真值轨迹(黑) + 站 1 检测云(蓝o) + 站 2 检测云经 R,t 变换后(红x),
%       坐标系 = Roomframe_plane 构造房间系(目标平面 z=0, 轨迹显示为水平);
%       nargin < 6 时退化为站 1 系原样
%   右: 逐配对帧残差(mm) + RMSE 线
% 返回 figure 句柄供 main_motion 导出存盘
%
% 输入: traj - Scene_trajectory; P_hist/T_hist - Det_to_points 输出 (Nf,n_st,3)
%       both (Nf,1) 双站同时检出掩码; fuse - Fuse_SVD 输出
%       R1 (3,3)/p1 (3,1) 可选, Roomframe_plane 输出, 站 1 系 -> 构造房间系

    P1  = reshape(P_hist(both, 1, :), [], 3);           % 站 1 检测云(站 1 系)
    P2t = reshape(P_hist(both, 2, :), [], 3) * fuse.R.' + fuse.t.';  % 站 2 -> 站 1
    T1  = reshape(T_hist(both, 1, :), [], 3);           % 真值轨迹(站 1 系)
    t_ms = traj.t(both) * 1e3;

    % ---- 站 1 系 -> 构造房间系(缺省 = 恒等, 即站 1 系原样) ----
    if nargin < 6, R1 = eye(3);  p1 = zeros(3,1); end
    P1r  = P1  * R1.' + p1.';
    P2tr = P2t * R1.' + p1.';
    T1r  = T1  * R1.' + p1.';

    fig = figure('Name', '双站点云融合', 'Color', 'w', 'Position', [80 80 1050 460]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('SVD 刚体配准(站 2 -> 站 1), 配对 %d 帧, RMSE %.1f mm', ...
          size(P1, 1), fuse.rmse*1e3));

    % ---- 左: 构造房间系点云叠加 ----
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    plot3(ax, T1r(:,1), T1r(:,2), T1r(:,3), 'k-', 'LineWidth', 1.4, 'DisplayName', '真值轨迹');
    plot3(ax, P1r(:,1), P1r(:,2), P1r(:,3), 'o', 'MarkerSize', 4, 'Color', [0.10 0.30 0.85], ...
          'DisplayName', '站 1 检测云');
    plot3(ax, P2tr(:,1), P2tr(:,2), P2tr(:,3), 'x', 'MarkerSize', 5, 'Color', [0.85 0.10 0.10], ...
          'LineWidth', 1.2, 'DisplayName', '站 2 检测云(R,t 变换后)');
    axis(ax, 'equal');
    xlabel(ax, 'x_w (m)'); ylabel(ax, 'y_w (m)'); zlabel(ax, 'z_w (m)');
    title(ax, '构造房间系点云配准叠加(目标平面 z = 0)');
    legend(ax, 'Location', 'best', 'FontSize', 8);
    view(ax, 40, 20);

    % ---- 右: 配对残差(全配对帧, 含 3σ 外点) ----
    if isfield(fuse, 'res_all'), res = fuse.res_all; else, res = fuse.res; end
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, t_ms, res * 1e3, 'b-', 'LineWidth', 1.1, 'DisplayName', '逐帧残差');
    yline(ax, fuse.rmse * 1e3, 'r--', 'LineWidth', 1.2, ...
          'DisplayName', sprintf('内点 RMSE %.1f mm', fuse.rmse*1e3));
    xlabel(ax, 't (ms)'); ylabel(ax, '残差 (mm)');
    title(ax, '配对残差(变换后站 2 点 vs 站 1 点)');
    legend(ax, 'Location', 'best', 'FontSize', 8);
end
