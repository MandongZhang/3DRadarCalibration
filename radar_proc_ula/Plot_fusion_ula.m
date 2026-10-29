function fig = Plot_fusion_ula(traj, T1, fuse, both)
% ULA 配准结果图: 左 = 拟合平面内真值轨迹与两站交点云(站1系基投影);
%                右 = 逐配对帧残差
% 输入: traj; T1 (Nf,3) 真值轨迹(站1真位姿系); fuse(Register_iterSVD 输出);
%       both (Nf,1) 配对掩码
    nu = fuse.plane1.nu;
    c0 = mean(fuse.Cs1, 1);
    e1 = [1;0;0] - nu(1)*nu;   e1 = e1/norm(e1);     % 面内基 1(去法向分量)
    e2 = cross(nu.', e1);                             % 面内基 2
    proj = @(P) [P*e1.', P*e2.'];
    T1b = T1(both, :);
    Tp = proj(T1b - ((T1b - c0)*nu.').*nu);           % 真值正交投影到拟合平面
    P1 = proj(fuse.Cs1);
    P2 = proj(fuse.Cs2);
    t_ms = traj.t(both)*1e3;
    fig = figure('Name', 'ULA 双站配准', 'Color', 'w', 'Position', [80 80 1050 460]);
    tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('ULA 平面约束迭代配准(站 2 -> 站 1), 配对 %d 帧, RMSE %.1f mm', ...
          size(fuse.Cs1, 1), fuse.rmse*1e3));
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    plot(ax, Tp(:,1), Tp(:,2), 'k-', 'LineWidth', 1.4, 'DisplayName', '真值轨迹(平面投影)');
    plot(ax, P1(:,1), P1(:,2), 'o', 'MarkerSize', 4, 'Color', [0.10 0.30 0.85], ...
         'DisplayName', '站 1 交点云');
    plot(ax, P2(:,1), P2(:,2), 'x', 'MarkerSize', 5, 'Color', [0.85 0.10 0.10], ...
         'LineWidth', 1.2, 'DisplayName', '站 2 交点云(配准后)');
    axis(ax, 'equal');
    xlabel(ax, '面内 e_1 (m)'); ylabel(ax, '面内 e_2 (m)');
    title(ax, '拟合平面内点云配准叠加');
    legend(ax, 'Location', 'best', 'FontSize', 8);
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plot(ax, t_ms, fuse.res*1e3, 'b-', 'LineWidth', 1.1, 'DisplayName', '逐帧残差');
    yline(ax, fuse.rmse*1e3, 'r--', 'LineWidth', 1.2, ...
          'DisplayName', sprintf('RMSE %.1f mm', fuse.rmse*1e3));
    xlabel(ax, 't (ms)'); ylabel(ax, '残差 (mm)');
    title(ax, '配对残差(变换后站 2 点 vs 站 1 点)');
    legend(ax, 'Location', 'best', 'FontSize', 8);
end
