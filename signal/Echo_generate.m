function echo = Echo_generate(cfg, scene, tgt, wf, noise_on)
% 信号级中频回波生成(第二块: 信号生成)——全仿真核心
%
% 输入: cfg; scene - Scene_radarPoses 输出; tgt - Scene_targets 输出;
%       wf - Waveform_derive 输出; noise_on - 可选, 默认 true, false 生成
%       无噪回波(自检用, 验证模型确定性部分)
% 信号级中频回波生成(第二块: 信号生成)——全仿真核心
%
% 输入: cfg; scene - Scene_radarPoses 输出; tgt - Scene_targets 输出;
%       wf - Waveform_derive 输出
% 输出: echo(i) - 每站一份的结构体数组
%   echo(i).data       (N_fast, ntx*N_chirps, nrx) 复双精度 IF 采样,
%                      噪底归一化为单位复功率(噪底 0 dB, 绝对量经 noise_dbm 换算)
%   echo(i).tx_of_slot (ntx*N_chirps,1)  slot->TX 映射 mod(s-1,ntx)+1
%   echo(i).virt       (3, ntx*nrx) 虚拟阵元房间系坐标, ch = (p-1)*nrx + q
%   echo(i).t_fast / .fc / .S / .lambda / .Fs / .N_chirps   波形透传
%   echo(i).noise_dbm  物理噪底参考(绝对 SNR 换算用)
%   echo(i).n_facet    (1,ntx) 每 TX 实际参与叠加的面元数(诊断)
%
% 模型(Echo_Simulation_Design.md, 双基近场精确, 逐面元 m / TX p / RX q):
%   tau_mpq = (||p_m - p_tx_p|| + ||p_rx_q - p_m||)/c        式 E3
%   s(t) = sum_m a_mpq exp( j2pi( fc*tau + S*tau*t - S*tau^2/2 ) )   式 E4
%   a_mpq = sqrt( Pt*Gtx*Grx*lambda^2*sigma_m / ((4pi)^3 Rtx^2 Rrx^2) / P_noise )
%           * exp(j*psi_m)                                    式 E9
%   sigma_m = (4pi/lambda^2) A_f^2 cos^2(theta_b)            式 E6, PO 平板
%   可见性: n·u_tx > 0 且 n·u_rx > 0 且方向图非零             式 E13
% 静止目标: 同 TX 各 chirp 确定部分相同, 信号合成一次复制到 chirp 维,
%           噪声逐 slot 独立(式 E5)
% 不建模: 直耦泄漏 / 站间互扰 / 多径 / chirp 边缘效应 / 接收链幅相不一致
% 调用方须先 rng(cfg.sim.seed)

    c = 299792458;
    if nargin < 5, noise_on = true; end
    N_fast = wf.N_fast;
    t = wf.t_fast;
    lambda = wf.lambda;

    % ---- 面元聚合(全部目标) ----
    M_tot = sum([tgt.M]);
    P   = zeros(3, M_tot);   Nrm = zeros(3, M_tot);
    Af  = zeros(M_tot, 1);   phf = zeros(M_tot, 1);
    m0 = 0;
    for k = 1:numel(tgt)
        idx = m0+1 : m0+tgt(k).M;
        P(:, idx)   = tgt(k).facet.pos.';
        Nrm(:, idx) = tgt(k).facet.nrm.';
        Af(idx)     = tgt(k).facet.area;
        phf(idx)    = tgt(k).facet.phase;
        m0 = m0 + tgt(k).M;
    end
    Nrm3 = reshape(Nrm, 3, M_tot, 1);

    % ---- 雷达方程常数(线性) ----
    Pt_lin      = 10^(cfg.radar(1).Pt_dbm/10) * 1e-3;      % W
    G0_lin      = 10^(cfg.radar(1).G0_dbi/10);
    Pnoise_lin  = 10^((wf.noise_dbm - 30)/10);             % W (dBm -> W 减 30)

    for i = 1:numel(scene.radar)
        rad  = scene.radar(i);
        ntx  = rad.ntx;   nrx = rad.nrx;
        n_ch = cfg.waveform.N_chirps;
        S_slots = ntx * n_ch;
        data = complex(zeros(N_fast, S_slots, nrx));
        n_facet = zeros(1, ntx);

        for p = 1:ntx
            % ---- TX 侧几何: (M,1) ----
            % uTx 统一定义为 面元->TX 方向: n·uTx > 0 = 法向朝 TX = 被照明
            % 方向图用其反方向(阵元指向面元), 两者勿混
            dTx = rad.tx_pos(:, p) - P;                    % (3,M)
            Rtx = vecnorm(dTx, 2, 1).';                    % (M,1)
            uTx = dTx ./ Rtx.';                            % (3,M) 面元->TX
            nTx = reshape(sum(Nrm .* uTx, 1), M_tot, 1);   % n·u_tx (M,1)
            Gtx = G0_lin * Antenna_gain(-uTx, rad.boresight, rad.fov, rad.g_p);

            % ---- RX 侧几何: (M,nrx) ----
            dRx = reshape(rad.rx_pos, 3, 1, nrx) - reshape(P, 3, M_tot, 1);
            Rrx = reshape(vecnorm(dRx, 2, 1), M_tot, nrx);
            uRx = dRx ./ reshape(Rrx, 1, M_tot, nrx);      % (3,M,nrx) 面元->RX
            nRx = reshape(sum(Nrm3 .* uRx, 1), M_tot, nrx);
            Grx = G0_lin * reshape(Antenna_gain(reshape(-uRx, 3, []), ...
                     rad.boresight, rad.fov, rad.g_p), M_tot, nrx);

            % ---- 可见性(式 E13)与 PO 平板 RCS(式 E6) ----
            vis = (nTx > 0 & Gtx > 0) & (nRx > 0 & Grx > 0);       % (M,nrx)
            uB  = uTx + uRx;                               % (3,M,nrx) 平分线方向
            nB  = reshape(vecnorm(uB, 2, 1), M_tot, nrx);
            cosb = max((nTx + nRx) ./ max(nB, eps), 0);    % clamp >= 0
            sigma = (4*pi/lambda^2) * (Af.^2) .* cosb.^2;  % (M,nrx)

            % ---- 双基雷达方程(式 E7)与复幅度(式 E9) ----
            Prx = Pt_lin .* (Gtx .* Grx) .* lambda^2 .* sigma ...
                  ./ ((4*pi)^3 .* (Rtx.^2) .* (Rrx.^2));   % (M,nrx) W
            W = sqrt(max(Prx, 0) / Pnoise_lin) .* vis .* exp(1i * phf);

            % ---- 双基延迟与相位(式 E3 E4) ----
            TAU  = (Rtx + Rrx) / c;                        % (M,nrx)
            beta = exp(1i * 2*pi * (wf.fc * TAU - wf.S/2 .* TAU.^2));
            fb   = wf.S .* TAU;                            % 拍频 (M,nrx)

            % ---- 快时间合成(面元分块广播, 控制临时量 ~8 MB) ----
            s0 = complex(zeros(N_fast, nrx));
            Mblk = 512;
            for m0b = 1:Mblk:M_tot
                idx = m0b : min(m0b + Mblk - 1, M_tot);
                Mb  = numel(idx);
                Wk  = W(idx, :) .* beta(idx, :);           % (Mb,nrx)
                Eb  = exp(1i*2*pi * reshape(fb(idx,:), Mb, 1, nrx) ...
                          .* reshape(t, 1, N_fast, 1));    % (Mb,N_fast,nrx)
                s0 = s0 + reshape(sum(Eb .* reshape(Wk, Mb, 1, nrx), 1), ...
                                  N_fast, nrx);
            end

            % ---- TDM 堆叠(式 E5): 信号复制, 噪声逐 slot 独立 ----
            cols = p : ntx : S_slots;
            s0_3d = reshape(s0, N_fast, 1, nrx);           % 先升维再复制
            if noise_on
                noise = (randn(N_fast, n_ch, nrx) + ...
                         1i*randn(N_fast, n_ch, nrx)) / sqrt(2);
                data(:, cols, :) = repmat(s0_3d, 1, n_ch, 1) + noise;
            else
                data(:, cols, :) = repmat(s0_3d, 1, n_ch, 1);
            end

            n_facet(p) = nnz(any(vis, 2));
        end

        % ---- 元数据 ----
        echo(i).data = data;
        echo(i).tx_of_slot = mod(0:S_slots-1, ntx)' + 1;
        echo(i).virt = zeros(3, ntx*nrx);
        for pp = 1:ntx
            for qq = 1:nrx
                if strcmp(cfg.array.type, 'ULA')
                    % ULA: 虚拟元沿站系 y, 连续半波长栅格(TX 长距 + RX 密集)
                    v_s = [(pp-1)*rad.nrx*rad.d_rx + (qq-1)*rad.d_rx; 0; 0];
                else
                    % L 型: TX 竖臂沿 z, RX 横臂沿 y
                    v_s = [0; 0; (pp-1)*rad.d_tx];
                    v_s(2) = qq*rad.d_rx;
                end
                echo(i).virt(:, (pp-1)*nrx + qq) = rad.p_true.' + rad.R_true * v_s;
            end
        end
        echo(i).t_fast = t;        echo(i).fc = wf.fc;
        echo(i).S = wf.S;          echo(i).lambda = lambda;
        echo(i).Fs = wf.Fs;        echo(i).N_chirps = n_ch;
        echo(i).d_tx = rad.d_tx;   echo(i).d_rx = rad.d_rx;   % 角度换算用
        echo(i).noise_dbm = wf.noise_dbm;
        echo(i).n_facet = n_facet;
    end
end
