function det = Process_cube_ula(cfg, scene, echo)
% ULA 处理链主入口: RangeFFT(ref 复用) + AngleFFT_ula
% det(i).points = [r (m) | u_y | 峰SNR (dB)]; det(i).spec1d 供绘图
    for i = 1:numel(echo)
        rad = scene.radar(i);
        [Rcube, r_ax, P_acc] = RangeFFT(echo(i));
        % 距离: 谱最大峰 + 抛物线内插(与 ref 同款)
        [~, k] = max(P_acc);
        dr = r_ax(2) - r_ax(1);
        a = P_acc(k-1);  b = P_acc(k);  c2 = P_acc(k+1);
        den = a - 2*b + c2;
        if den ~= 0, dk = 0.5*(a-c2)/den; else, dk = 0; end
        dk = max(min(dk,1),-1);
        r_est = r_ax(k) + dk*dr;
        % 1D 测角: 原始观测量 u_y = sin(az)cos(el)
        [uy, S1db, ax] = AngleFFT_ula(Rcube, k, rad.ntx, echo(i).lambda, ...
                                      echo(i).d_rx, cfg.proc.nfft_az);
        det(i).points = [r_est, uy, 10*log10(P_acc(k)/median(P_acc))];
        det(i).r_ax = r_ax;   det(i).P_acc = P_acc;
        det(i).spec1d = struct('S', S1db, 'x', ax);
    end
end
