
% 3D 双雷达轨迹遍历仿真 - 运动目标入口(Z/C/M 静态点遍历)
%
% 思路: 轨迹看成多个静态点的序列。外层循环逐帧把目标平移到轨迹点, 每帧内部
% 与静态流程完全一致(Echo_generate -> Process_cube -> Truth_project),
% 循环内部不画图; 循环结束后按 cfg.sim.plot 开关统一出图并存盘。
% 停-走冻结近似与多普勒可忽略论证: Echo_Simulation_Design.md 运动目标节。
% 运行: MATLAB 中运行 main_motion, 或命令行 matlab -batch "main_motion"

    fprintf('==============================================\n');
    fprintf('   3D 双雷达轨迹遍历仿真(运动目标)\n');
    fprintf('==============================================\n\n');

    addpath('scene', 'signal', 'radar_proc');
    cfg   = config();
    if ~strcmp(cfg.array.type, 'L')
        error('main_motion 为 ref(L 型)入口: 请将 cfg.array.type 设为 ''L'', 或改用 main_motion_ula');
    end
    rng(cfg.sim.seed);

    scene = Scene_radarPoses(cfg);
    tgt0  = Scene_targets(cfg);
    wf    = Waveform_derive(cfg);
    traj  = Scene_trajectory(cfg, wf);

    Nf    = size(traj.pos, 1);
    n_st  = numel(scene.radar);
    T_frame = wf.T_frame;

    fprintf('【场景】房间 %.1f x %.1f x %.1f m, 雷达 %d 站\n', ...
            cfg.room.length, cfg.room.width, cfg.room.height, n_st);
    fprintf('【轨迹】%s 型, %.1f m 高度, 速度 %.1f m/s, %d 帧\n', ...
            traj.type, cfg.target.targets(1).pos(3), cfg.target.traj.speed, Nf);
    fprintf('  帧周期 %.2f ms = 采集 %.2f + 空闲 %.2f, 实际帧率 %.1f Hz, 总时长 %.2f s\n', ...
            wf.T_frame*1e3, (wf.T_frame - wf.T_idle)*1e3, wf.T_idle*1e3, ...
            wf.frame_rate, traj.t(end));
    fprintf('  起点 [%.2f %.2f %.2f] -> 终点 [%.2f %.2f %.2f]\n\n', ...
            traj.pos(1,:), traj.pos(end,:));

    % ---- 轨迹遍历主循环: 每帧 = 一个静态点(循环内不画图) ----
    pt_hist    = nan(Nf, n_st, 4);        % [r | az | el | 峰SNR] 逐帧逐站
    truth_hist = nan(Nf, n_st, 3);        % [r | az | el] 目标中心真值投影
    detected   = false(Nf, n_st);         % 峰SNR >= snr_th
    det_last   = [];   truth_last = [];   % 末帧留作出图

    tic;
    for n = 1:Nf
        tgt_n  = Move_target(tgt0, traj.pos(n,:) - tgt0(1).pos);
        echo   = Echo_generate(cfg, scene, tgt_n, wf);
        det    = Process_cube(cfg, scene, echo);
        truth  = Truth_project(scene, tgt_n);
        for i = 1:n_st
            pt_hist(n, i, :)    = det(i).points;
            truth_hist(n, i, :) = truth(i).truth(1, :);
            detected(n, i)      = det(i).points(4) >= cfg.motion.snr_th;
        end
        det_last = det;   truth_last = truth;

        if mod(n, 10) == 0 || n == Nf
            fprintf('  帧 %3d/%d  t = %6.1f ms  pos = [%.2f %.2f %.2f]  SNR(dB):', ...
                    n, Nf, traj.t(n)*1e3, traj.pos(n,:));
            for i = 1:n_st
                fprintf('  站%d %5.1f%s', i, det(i).points(4), ternary(detected(n,i), '✓', '✗'));
            end
            fprintf('\n');
        end
    end
    t_elapsed = toc;
    fprintf('  遍历完成: %d 帧, 耗时 %.1f s (%.0f ms/帧)\n\n', Nf, t_elapsed, t_elapsed/Nf*1e3);

    % ---- 汇总统计(检出帧上的误差) ----
    fprintf('【汇总】检测判定阈值 %d dB\n', cfg.motion.snr_th);
    for i = 1:n_st
        d = detected(:, i);
        err   = squeeze(pt_hist(:, i, 1:3)) - squeeze(truth_hist(:, i, 1:3));   % (Nf,3)
        snr_i = squeeze(pt_hist(:, i, 4));
        fprintf('  雷达 %d: 检出 %d/%d 帧 (%.0f%%)', i, nnz(d), Nf, 100*nnz(d)/Nf);
        if any(d)
            fprintf(', 检出帧误差均值 |Δr| = %.3f m, |Δaz| = %.2f deg, |Δel| = %.2f deg, 峰SNR中位 %.1f dB', ...
                    mean(abs(err(d,1))), mean(abs(err(d,2))), mean(abs(err(d,3))), ...
                    median(snr_i(d)));
        end
        fprintf('\n');
    end

    %% 点云累积 + 双站 SVD 融合(自校准核心: 估计站 2 -> 站 1 相对外参)
    % 逐帧检测 (r,az,el) -> 各站站系直角坐标, 帧序即配对关系; 双站同时检出的
    % 帧进入配准。真值关系: C_s1 = (R1'R2) C_s2 + R1'(p2 - p1)
    P_hist = Det_to_points(pt_hist);           % (Nf,n_st,3) 检测点云(各站系)
    T_hist = Det_to_points(truth_hist);        % 真值点云(绘图对照)
    both   = all(detected, 2);                 % 双站同时检出 = 配对帧
    P1 = reshape(P_hist(both, 1, :), [], 3);   % 站 1 系配对点云
    P2 = reshape(P_hist(both, 2, :), [], 3);   % 站 2 系配对点云
    fuse1 = Fuse_SVD(P1, P2);                  % 第一遍: 全配对解
    inl   = fuse1.res < 3 * fuse1.rmse;        % 3σ 剔除大残差配对帧
    fuse  = Fuse_SVD(P1(inl, :), P2(inl, :));  % 第二遍: 内点重解
    fuse.res_all = fuse1.res;                  % 全配对帧残差(绘图用, 含外点)
    fuse.n_inlier = nnz(inl);                  % 内点数

    R_rel = scene.radar(1).R_true.' * scene.radar(2).R_true;
    t_rel = scene.radar(1).R_true.' * (scene.radar(2).p_true - scene.radar(1).p_true).';
    eR = acosd(max(min((trace(R_rel.' * fuse.R) - 1)/2, 1), -1));
    et = norm(fuse.t - t_rel);
    fprintf('【融合】SVD 刚体配准(站 2 -> 站 1), 配对 %d 帧, 3σ 内点 %d (%.0f%%)\n', ...
            nnz(both), fuse.n_inlier, 100*fuse.n_inlier/nnz(both));
    fprintf('  旋转误差 %.4f deg, 平移误差 %.2f mm, 内点 RMSE %.2f mm\n', ...
            eR, et*1e3, fuse.rmse*1e3);
    fprintf('  t_hat = [%+.4f %+.4f %+.4f] m, t_true = [%+.4f %+.4f %+.4f] m\n', ...
            fuse.t, t_rel);

    % ---- 构造房间系(平面拟合, 仅用于显示): z=0 取目标平面, 原点在站 1 正下,
    %      x = 站 1 波束水平投影; 雷达离面高度 d 由拟合直接可观 ----
    [R1c, p1c, pf] = Roomframe_plane(P1(inl, :));
    fprintf('  构造房间系: 雷达离目标平面 %.3f m, 平面拟合 RMS %.2f mm\n', ...
            pf.d, pf.rms*1e3);

    %% 出图(cfg.sim.plot 开关; 循环已结束, 此处统一画)
    if cfg.sim.plot
        figT = Plot_tracking(cfg, scene, traj, pt_hist, truth_hist, detected);
        figP = Plot_processing(det_last, truth_last, cfg);
        figF = Plot_fusion(traj, P_hist, T_hist, both, fuse, R1c, p1c);
    end

    %% 存结果: 轨迹与逐帧检测记录 + 点云与融合结果(不存逐帧 cube, 省内存)
    if cfg.sim.save
        if ~exist('results', 'dir'), mkdir('results'); end
        tag = char(datetime('now', 'Format', 'yyyy-MM-dd'));
        if cfg.sim.plot
            exportgraphics(figT, fullfile('results', [tag '_tracking.png']), 'Resolution', 150);
            exportgraphics(figP, fullfile('results', [tag '_motion_proc.png']), 'Resolution', 150);
            exportgraphics(figF, fullfile('results', [tag '_fusion.png']), 'Resolution', 150);
        end
        save(fullfile('results', [tag '_motion.mat']), 'cfg', 'traj', 'T_frame', ...
             'pt_hist', 'truth_hist', 'detected', ...
             'P_hist', 'T_hist', 'both', 'fuse', 'R_rel', 't_rel');
        if cfg.sim.plot
            fprintf('\n结果已存 results/ (%s_tracking.png, %s_motion_proc.png, %s_fusion.png, %s_motion.mat)\n', tag, tag, tag, tag);
        else
            fprintf('\n结果已存 results/ (%s_motion.mat, plot 关闭未出图)\n', tag);
        end
    end

% ================================================================
function s = ternary(cond, a, b)
% 标量条件取值(打印标记用)
    if cond, s = a; else, s = b; end
end
