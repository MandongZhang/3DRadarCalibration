function fig = Plot_processing(det, truth, cfg)
% 处理链结果图(tiledlayout 2x2): 每站一行 = 距离谱 + 2D 角度谱
% 角度谱显示为角度空间: 方向余弦谱(FFT 原生, 线性 u 刻度)经 interp2 重采样
% 到均匀 (az, el) 网格, 横轴方位角 az、纵轴俯仰角 el (deg), 星标 = 检测角
% 显示范围 = +-cfg.radar(i).fov(半锥角, 与天线方向图支撑域一致); 锥外留白,
% 黑虚线 = fov 边界(锥面投影)
% 不传 cfg 时退化为默认 60 deg
% 返回 figure 句柄供 main 导出存盘

    n_st = numel(det);
    fig = figure('Name', '处理链结果', 'Color', 'w', 'Position', [60 60 980 620]);
    tl = tiledlayout(fig, n_st, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, 'FFT 处理链: 距离谱与角度谱');

    for i = 1:n_st
        % ---- 距离谱 ----
        ax1 = nexttile((i-1)*2 + 1);
        hold(ax1, 'on'); grid(ax1, 'on');
        plot(ax1, det(i).r_ax, 10*log10(det(i).P_acc), 'b-', 'LineWidth', 1.2, ...
             'DisplayName', '非相干积累距离谱');
        hold(ax1, 'off');
        xline(ax1, det(i).points(1), 'r--', 'LineWidth', 1.3, ...
              'DisplayName', sprintf('检测 r = %.2f m', det(i).points(1)));
        ylim(ax1, 'auto');
        xlabel(ax1, '距离 (m)'); ylabel(ax1, '功率 (dB)');
        title(ax1, sprintf('雷达 %d 距离谱', i));
        legend(ax1, 'Location', 'northeast', 'FontSize', 8);

        % ---- 2D 角度谱(方向余弦谱 -> 角度空间重采样) ----
        ax2 = nexttile((i-1)*2 + 2);
        sp = det(i).spec2d;
        if nargin < 3, fov_i = 60; else, fov_i = cfg.radar(i).fov; end
        az_lim = fov_i;  el_lim = fov_i;           % 显示范围 (deg) = 半锥角 fov
        az_g = linspace(-az_lim, az_lim, 361);
        el_g = linspace(-el_lim, el_lim, 361);
        [AZ, EL] = meshgrid(az_g, el_g);           % (n_el, n_az)
        Uy = sind(AZ) .* cosd(EL);                 % 横轴方向余弦 sin(az)cos(el)
        Uz = sind(EL);                             % 纵轴方向余弦 sin(el)
        S_ang = interp2(sp.x, sp.y, sp.S', Uy, Uz, 'linear');   % (n_el, n_az)
        S_ang((Uy.^2 + Uz.^2) > sind(fov_i)^2) = NaN;   % 锥外增益为零, 留白
        imagesc(ax2, az_g, el_g, S_ang);
        hold(ax2, 'on');
        phi = linspace(0, 2*pi, 361);                   % 锥角边界(方向余弦圆的 az/el 投影)
        sf  = sind(fov_i);
        el_b = asind(max(min(sf .* sin(phi), 1), -1));
        arg  = sf .* cos(phi) ./ sqrt(max(1 - (sf .* sin(phi)).^2, 0));
        az_b = asind(max(min(arg, 1), -1));   % 钳制±1: fov=90 时 0/0->NaN 由 plot 跳过, 防浮点越界出复数
        plot(ax2, az_b, el_b, 'k--', 'LineWidth', 1.2);
        pt = det(i).points;
        plot(ax2, pt(2), pt(3), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
        hold(ax2, 'off');
        xlim(ax2, [-az_lim, az_lim]);  ylim(ax2, [-el_lim, el_lim]);   % 钉死显示窗 = ±fov
        set(ax2, 'YDir', 'normal');
        xlabel(ax2, '方位角 az (deg)'); ylabel(ax2, '俯仰角 el (deg)');
        title(ax2, sprintf('归一化角度谱(黑虚线=fov 边界), 估计 az = %.1f°, el = %.1f°', pt(2), pt(3)));
        colormap(ax2, 'jet'); colorbar(ax2);
    end
end
