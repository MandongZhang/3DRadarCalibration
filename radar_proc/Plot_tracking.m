function fig = Plot_tracking(cfg, scene, traj, pt_hist, truth_hist, detected)
% 轨迹遍历结果总图(tiledlayout 4x2, main_motion 循环结束后统一调用)
%   左上(跨2行): 3D 房间 + 真值轨迹 + 双站检测点(峰SNR 着色, 未检出灰叉)
%   右上(跨2行): 峰SNR 时间线 + 判定阈值
%   下两行:     r / az / el 估计 vs 真值 | 检测时间线
% 返回 figure 句柄供 main_motion 导出存盘
%
% 输入: cfg; scene - Scene_radarPoses; traj - Scene_trajectory
%       pt_hist    (Nf, n_st, 4) [r | az | el | 峰SNR(dB)] 逐帧逐站
%       truth_hist (Nf, n_st, 3) [r | az | el] 真值投影
%       detected   (Nf, n_st) logical

    Nf   = size(pt_hist, 1);
    n_st = size(pt_hist, 2);
    t_ms = traj.t * 1e3;
    snr  = squeeze(pt_hist(:, :, 4));                       % (Nf, n_st)
    cols = {'r', 'b'};
    T_frame = cfg.radar(1).ntx * cfg.waveform.N_chirps * cfg.waveform.Tc;

    fig = figure('Name', '轨迹遍历结果', 'Color', 'w', 'Position', [50 40 1220 880]);
    tl = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('轨迹遍历仿真: %s 型轨迹, %d 帧, 帧周期 %.1f ms', ...
          traj.type, Nf, T_frame*1e3));

    % ---- (1-2,1) 3D: 房间 + 轨迹 + 检测点 ----
    ax = nexttile([2 1]); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    L = cfg.room.length;  W = cfg.room.width;  H = cfg.room.height;
    plot3(ax, [0 L L 0 0], [0 0 W W 0], [0 0 0 0 0], 'k-', 'LineWidth', 1, ...
          'HandleVisibility', 'off');
    for cx = [0 L]
        for cy = [0 W]
            plot3(ax, [cx cx], [cy cy], [0 H], 'k--', 'HandleVisibility', 'off');
        end
    end
    plot3(ax, [0 L L 0 0], [0 0 W W 0], [H H H H H], 'k:', 'LineWidth', 0.5, ...
          'HandleVisibility', 'off');

    plot3(ax, traj.pos(:,1), traj.pos(:,2), traj.pos(:,3), 'k-', 'LineWidth', 1.6, ...
          'DisplayName', '真值轨迹');
    plot3(ax, traj.pos(1,1), traj.pos(1,2), traj.pos(1,3), 'k>', 'MarkerFaceColor', 'k', ...
          'MarkerSize', 8, 'DisplayName', '起点');
    plot3(ax, traj.pos(end,1), traj.pos(end,2), traj.pos(end,3), 'ks', ...
          'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', '终点');

    for i = 1:n_st                                          % 雷达 + boresight
        p = scene.radar(i).p_true;  b = scene.radar(i).boresight;
        scatter3(ax, p(1), p(2), p(3), 90, cols{i}, 'p', 'filled', ...
                 'DisplayName', sprintf('雷达 %d', i));
        quiver3(ax, p(1), p(2), p(3), b(1)*2, b(2)*2, b(3)*2, 0, ...
                'Color', cols{i}, 'LineWidth', 1.4, 'HandleVisibility', 'off');
    end

    for i = 1:n_st                                          % 检测点(站系球坐标 -> 房间系)
        est = squeeze(pt_hist(:, i, 1:3));                  % (Nf,3) [r az el]
        Cw  = zeros(Nf, 3);
        for n = 1:Nf
            Cw(n, :) = spher2room(est(n, :), scene.radar(i));
        end
        d = detected(:, i);
        if any(d)
            scatter3(ax, Cw(d,1), Cw(d,2), Cw(d,3), 26, snr(d,i), 'filled', ...
                     'DisplayName', sprintf('站 %d 检出', i));
        end
        if any(~d)
            scatter3(ax, Cw(~d,1), Cw(~d,2), Cw(~d,3), 30, [0.55 0.55 0.55], 'x', ...
                     'DisplayName', '未检出');
        end
    end

    xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)'); zlabel(ax, 'z (m)');
    axis(ax, 'equal'); xlim(ax, [-0.5 L+0.5]); ylim(ax, [-0.5 W+0.5]); zlim(ax, [0 H+0.5]);
    view(ax, 38, 22);
    colormap(ax, 'jet');
    cb = colorbar(ax);  cb.Label.String = '峰SNR (dB)';  cb.Label.FontSize = 8;
    title(ax, '3D 轨迹与双站检测点(SNR 着色)');
    legend(ax, 'Location', 'northeastoutside', 'FontSize', 8);

    % ---- (1-2,2) SNR 时间线 ----
    ax = nexttile([2 1]); hold(ax, 'on'); grid(ax, 'on');
    for i = 1:n_st
        plot(ax, t_ms, snr(:, i), '-', 'Color', cols{i}, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('站 %d 峰SNR', i));
        plot(ax, t_ms(detected(:,i)), snr(detected(:,i), i), '.', 'Color', cols{i}, ...
             'MarkerSize', 12, 'HandleVisibility', 'off');
    end
    yline(ax, cfg.motion.snr_th, 'k--', 'LineWidth', 1.1, ...
          'DisplayName', sprintf('阈值 %d dB', cfg.motion.snr_th));
    xlabel(ax, 't (ms)'); ylabel(ax, '峰SNR (dB)');
    title(ax, '检测 SNR 时间线');
    legend(ax, 'Location', 'best', 'FontSize', 8);

    % ---- (3,1)(3,2)(4,1) r / az / el: 估计 vs 真值 ----
    qty = {'r', 'az', 'el'};
    ylab = {'r (m)', 'az (deg)', 'el (deg)'};
    for q = 1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        for i = 1:n_st
            est_q   = squeeze(pt_hist(:, i, q));
            truth_q = squeeze(truth_hist(:, i, q));
            plot(ax, t_ms, est_q, '-',  'Color', cols{i}, 'LineWidth', 1.2, ...
                 'DisplayName', sprintf('站 %d 估计', i));
            plot(ax, t_ms, truth_q, '--', 'Color', cols{i}, 'LineWidth', 1.0, ...
                 'DisplayName', sprintf('站 %d 真值', i));
        end
        xlabel(ax, 't (ms)'); ylabel(ax, ylab{q});
        title(ax, sprintf('%s 估计 vs 真值', qty{q}));
        if q == 1, legend(ax, 'Location', 'best', 'FontSize', 7, 'NumColumns', 2); end
    end

    % ---- (4,2) 检测时间线(白 = 未检出, 绿 = 检出) ----
    ax = nexttile;
    imagesc(ax, t_ms, 1:n_st, double(detected.'));
    set(ax, 'YDir', 'normal');
    yticks(ax, 1:n_st); yticklabels(ax, compose('站 %d', 1:n_st));
    colormap(ax, [1 1 1; 0.2 0.6 0.2]);
    xlabel(ax, 't (ms)');
    title(ax, '检测时间线(绿 = 检出)');
end

% ================================================================
function Cw = spher2room(pt, rad)
% 站系球坐标 (r, az, el) -> 房间系直角坐标(与 Truth_project 投影互逆)
    r = pt(1);  az = pt(2);  el = pt(3);
    Cs = r * [cosd(el)*cosd(az); cosd(el)*sind(az); sind(el)];
    Cw = rad.R_true * Cs + rad.p_true.';
end
