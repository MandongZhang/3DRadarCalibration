% 3D 双雷达 ULA 轨迹遍历仿真与平面约束迭代配准(新算法入口)
%
% 传感: ULA 雷达(TX/RX 共线横臂), 每站每帧只测 r 与 u_y = sin(az)cos(el);
% 场景/轨迹/真值外参与 ref(main_motion) 完全一致, 6D 外参为估计目标。
% 配准: Register_iterSVD - 平面约束圆弧求交 + Fuse_SVD 定点迭代,
%       初值为安装先验; 水平面假设(z = h, 物理事实)为算法前提。
% 对照: main_motion(ref) 的 <date>_motion.mat; Compare_ref_ula 并排比较。
% 运行: matlab -batch "main_motion_ula"

    fprintf('==============================================\n');
    fprintf('   3D 双雷达 ULA 轨迹遍历(平面约束迭代配准)\n');
    fprintf('==============================================\n\n');

    addpath('scene', 'signal', 'radar_proc', 'radar_proc_ula');
    cfg = config();
    cfg.array.type = 'ULA';                 % 本入口固定 ULA 阵列
    rng(cfg.sim.seed);

    scene = Scene_radarPoses(cfg);
    tgt0  = Scene_targets(cfg);
    wf    = Waveform_derive(cfg);
    traj  = Scene_trajectory(cfg, wf);

    Nf   = size(traj.pos, 1);
    n_st = numel(scene.radar);
    h    = cfg.target.targets(1).pos(3);    % 轨迹平面高度(算法前提)

    fprintf('【场景】房间 %.1f x %.1f x %.1f m, 雷达 %d 站 (阵列: %s)\n', ...
            cfg.room.length, cfg.room.width, cfg.room.height, n_st, cfg.array.type);
    fprintf('【轨迹】%s 型, 高度 %.2f m, 速度 %.1f m/s, %d 帧\n\n', ...
            traj.type, h, cfg.target.traj.speed, Nf);

    pt_hist    = nan(Nf, n_st, 3);          % [r | u_y | 峰SNR]
    truth_hist = nan(Nf, n_st, 2);          % [r | u_y] 真值
    detected   = false(Nf, n_st);

    tic;
    for n = 1:Nf
        tgt_n = Move_target(tgt0, traj.pos(n,:) - tgt0(1).pos);
        echo  = Echo_generate(cfg, scene, tgt_n, wf);
        det   = Process_cube_ula(cfg, scene, echo);
        for i = 1:n_st
            pt_hist(n, i, :) = det(i).points;
            Cs = scene.radar(i).R_true.' * (tgt_n(1).pos - scene.radar(i).p_true).';
            truth_hist(n, i, :) = [norm(Cs), Cs(2)/norm(Cs)];
            detected(n, i) = det(i).points(3) >= cfg.motion.snr_th;
        end
        if mod(n, 10) == 0 || n == Nf
            fprintf('  帧 %3d/%d  t = %6.1f ms  pos = [%.2f %.2f %.2f]  SNR(dB):', ...
                    n, Nf, traj.t(n)*1e3, traj.pos(n,:));
            for i = 1:n_st
                fprintf('  站%d %5.1f%s', i, det(i).points(3), ternary(detected(n,i), '✓', '✗'));
            end
            fprintf('\n');
        end
    end
    t_elapsed = toc;
    fprintf('  遍历完成: %d 帧, 耗时 %.1f s\n\n', Nf, t_elapsed);

    % ---- ULA 平面约束迭代配准(6D 外参估计) ----
    both = all(detected, 2);
    D1 = squeeze(pt_hist(both, 1, 1:2));
    D2 = squeeze(pt_hist(both, 2, 1:2));
    R1p = scene.radar(1).R_prior;  p1 = scene.radar(1).p_true;
    R2p = scene.radar(2).R_prior;  p2 = scene.radar(2).p_true;
    fuse = Register_iterSVD(D1, D2, R1p, p1, R2p, p2, h);
    R_rel = scene.radar(1).R_true.' * scene.radar(2).R_true;
    t_rel = scene.radar(1).R_true.' * (scene.radar(2).p_true - scene.radar(1).p_true).';
    eR = acosd(max(min((trace(R_rel.' * fuse.R) - 1)/2, 1), -1));
    et = norm(fuse.t - t_rel);
    fprintf('【融合】ULA 平面约束迭代配准(站 2 -> 站 1), 配对 %d 帧\n', nnz(both));
    fprintf('  旋转误差 %.4f deg, 平移误差 %.2f mm, 内点 RMSE %.2f mm, 迭代 %d 轮(收敛=%d)\n', ...
            eR, et*1e3, fuse.rmse*1e3, fuse.n_iter, fuse.conv);
    fprintf('  t_hat = [%+.4f %+.4f %+.4f] m, t_true = [%+.4f %+.4f %+.4f] m\n', ...
            fuse.t, t_rel);

    % ---- 真值轨迹(站1真位姿系)与出图/存盘 ----
    T1 = nan(Nf, 3);
    for n = 1:Nf
        T1(n,:) = (scene.radar(1).R_true.' * (traj.pos(n,:) - scene.radar(1).p_true).').';
    end
    if cfg.sim.plot
        figF = Plot_fusion_ula(traj, T1, fuse, both);
    end
    if cfg.sim.save
        if ~exist('results', 'dir'), mkdir('results'); end
        tag = char(datetime('now', 'Format', 'yyyy-MM-dd'));
        if cfg.sim.plot
            exportgraphics(figF, fullfile('results', [tag '_fusion_ula.png']), 'Resolution', 150);
        end
        save(fullfile('results', [tag '_motion_ula.mat']), 'cfg', 'traj', 'pt_hist', ...
             'truth_hist', 'detected', 'fuse', 'R_rel', 't_rel', 'T1');
        fprintf('\n结果已存 results/ (%s_motion_ula.mat)\n', tag);
    end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end
