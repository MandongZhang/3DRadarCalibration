
% 3D 双雷达自校准仿真 - 一键运行入口(信号级回波版本)
%
% 流程: config -> 场景 -> 目标面元 -> 波形 -> 每站回波 cube -> 自检 -> 出图 -> 存结果
% 运行: MATLAB 中运行 main, 或命令行 matlab -batch "main"
% 设计文档: Echo_Simulation_Design.md(式 E1-E13) 

    fprintf('==============================================\n');
    fprintf('   3D 双雷达自校准仿真(信号级回波)\n');
    fprintf('==============================================\n\n');

    addpath('scene', 'signal', 'radar_proc');
    cfg = config();
    rng(cfg.sim.seed);

    %% 第一块: 场景布局
    scene = Scene_radarPoses(cfg);
    tgt   = Scene_targets(cfg);
    wf    = Waveform_derive(cfg);

    % Print sence configration
    fprintf('【场景】房间 %.1f x %.1f x %.1f m\n', ...
            cfg.room.length, cfg.room.width, cfg.room.height);
    for i = 1:numel(scene.radar)
        fprintf('  雷达 %d: p_true = [%.2f %.2f %.2f], 阵列 %dTX x %dRX (L 型)\n', ...
                i, scene.radar(i).p_true, scene.radar(i).ntx, scene.radar(i).nrx);
    end
    fprintf('  目标: %d 个方体, 面元总数 %d (facet %.1f cm, %s 模式)\n', ...
            numel(tgt), sum([tgt.M]), cfg.target.facet_size*100, cfg.target.phase_mode);
    for k = 1:numel(tgt)
        fprintf('    方体 %d: pos = [%.2f %.2f %.2f], size = [%.2f %.2f %.2f]\n', ...
                k, tgt(k).pos, tgt(k).size);
    end
    fprintf('【波形】fc = %.1f GHz, B = %.1f GHz, S = %.1e Hz/s\n', ...
            wf.fc/1e9, wf.B/1e9, wf.S);
    fprintf('  N_fast = %d, N_chirps = %d/TX, 距离分辨率 %.1f cm, 噪底 %.1f dBm\n\n', ...
            wf.N_fast, wf.N_chirps, 100 * 3e8 / (2 * wf.B), wf.noise_dbm);

    %% 第二块: 每站独立生成回波 cube

    echo = Echo_generate(cfg, scene, tgt, wf);
    
    % Print Radar echo simulation
    for i = 1:numel(echo)
        sz = numel(echo(i).data) * 8 / 1e6;
        fprintf('【回波】雷达 %d cube: (%d, %d, %d) 复双精度, %.1f MB, 噪底 %.1f dBm\n', ...
                i, size(echo(i).data), sz, echo(i).noise_dbm);
        fprintf('  每 TX 参与面元数: [%s]\n', num2str(echo(i).n_facet));
    end
    fprintf('\n');

    %% 第三块: FFT 处理链 - 距离与角度
    det   = Process_cube(cfg, scene, echo);
    truth = Truth_project(scene, tgt);
    
    fprintf('【处理】FFT 处理链估计 vs 真值(站系视角)\n');
    for i = 1:numel(det)
        pt = det(i).points;
        tr = truth(i).truth;
        fprintf('  雷达 %d: 检测 r = %.3f m, az = %+.2f deg, el = %+.2f deg\n', ...
                i, pt(1), pt(2), pt(3));
        for k = 1:size(tr, 1)
            fprintf('    真值 %d: r = %.3f m, az = %+.2f deg, el = %+.2f deg, 误差(r/az/el) = %.3f m / %+.2f / %+.2f\n', ...
                    k, tr(k,1), tr(k,2), tr(k,3), ...
                    pt(1)-tr(k,1), pt(2)-tr(k,2), pt(3)-tr(k,3));
        end
    end
    fprintf('\n');

    %% 自检(六项)
    ok = Check_echo(cfg, scene, tgt, wf);

    %% 出图(cfg.sim.plot 开关: false 只算不画)
    if cfg.sim.plot
        fig1 = Plot_scene(cfg, scene, tgt);
        fig2 = Plot_antenna(cfg);
        fig3 = Plot_processing(det, truth, cfg);
    end

    %% 存结果: 每站 radar cube + 元数据 + 诊断图 + 处理结果
    if cfg.sim.save
        if ~exist('results', 'dir'), mkdir('results'); end
        tag = char(datetime('now', 'Format', 'yyyy-MM-dd'));
        if cfg.sim.plot
            exportgraphics(fig1, fullfile('results', [tag '_scene3d.png']), 'Resolution', 150);
            exportgraphics(fig2, fullfile('results', [tag '_antenna.png']), 'Resolution', 150);
            exportgraphics(fig3, fullfile('results', [tag '_radar_proc.png']), 'Resolution', 150);
        end
        save(fullfile('results', [tag '_echo.mat']), 'cfg', 'scene', 'tgt', 'wf', 'echo', ...
             'det', 'truth');
        if cfg.sim.plot
            fprintf('结果已存 results/ (%s_*.png, %s_echo.mat)\n', tag, tag);
        else
            fprintf('结果已存 results/ (%s_echo.mat, plot 关闭未出图)\n', tag);
        end
    end

    if ~ok
        error('存在未通过的自检项');
    end

