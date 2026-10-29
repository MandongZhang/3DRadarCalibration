function det = Process_cube(cfg, scene, echo)
% 信号处理链主入口(纯 FFT): 每站 echo cube -> 距离 + 角度
%
% 输入: cfg; scene - Scene_radarPoses; echo - Echo_generate 输出 (1xN 站)
% 输出: det(i) 结构体数组, 每站一份
%   det(i).points (1,4) [r (m) | az (deg) | el (deg) | 峰SNR (dB)]
%   det(i).r_ax / .P_acc  距离轴与非相干积累谱(绘图)
%   det(i).spec2d         fftshift 2D 角度谱(绘图)
%
% 链路: RangeFFT 取谱最大峰 -> 距离(抛物线内插); AngleFFT(补零 cfg.proc)
% -> (az, el)。当前取单目标(谱最大峰); 多目标时扩展为多峰搜索即可。

    for i = 1:numel(echo)
        rad = scene.radar(i);
        [Rcube, r_ax, P_acc] = RangeFFT(echo(i));

        % ---- 距离: FFT 谱最大峰 + 抛物线内插 ----
        [~, k] = max(P_acc);
        dr = r_ax(2) - r_ax(1);
        a = P_acc(k-1);  b = P_acc(k);  c = P_acc(k+1);
        den = a - 2*b + c;
        if den ~= 0, dk = 0.5*(a-c)/den; else, dk = 0; end
        dk = max(min(dk, 1), -1);
        r_est = r_ax(k) + dk * dr;

        % ---- 角度: 该 bin 的 2D 角度 FFT(补零) ----
        [az, el, S2db, ax_x, ax_y] = AngleFFT(Rcube, k, rad.ntx, ...
            echo(i).lambda, echo(i).d_tx, echo(i).d_rx, ...
            cfg.proc.nfft_az, cfg.proc.nfft_el);

        det(i).points = [r_est, az, el, 10*log10(P_acc(k)/median(P_acc))];
        det(i).r_ax = r_ax;
        det(i).P_acc = P_acc;
        det(i).spec2d = struct('S', S2db, 'x', ax_x, 'y', ax_y);
    end
end
