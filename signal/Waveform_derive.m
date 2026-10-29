function wf = Waveform_derive(cfg)
% FMCW 波形推导与合法性校验(第二块: 信号生成)
%
% 输入: cfg - config.m 输出
% 输出: wf - 波形结构体
%   wf.S         斜率 B/Tc (Hz/s)
%   wf.Ts        快时间采样间隔 (s)
%   wf.t_fast    (N_fast,1) 快时间轴 (s)
%   wf.N_fast    每 chirp 采样数
%   wf.lambda    波长 (m)
%   wf.fc / .B / .Tc / .Fs / .N_chirps   透传
%   wf.noise_dbm 物理噪底(式 E8, kT0 = -174 dBm/Hz @ 290K, 带宽取 Fs)
%   wf.T_idle    帧间空闲 (s), 由 cfg.waveform.frame_rate 推导, 0 = 背靠背
%   wf.T_frame   帧周期 (s) = 采集时长 ntx*N_chirps*Tc + T_idle
%   wf.frame_rate 实际帧率 (Hz) = 1/T_frame
%
% 校验(失败即 error, 式 E2):
%   最大拍频 S*(2*max_range + 阵列口径) < Fs/2   (双基总路径)
%   Tc*Fs 须接近整数(N_fast)

    c = 299792458;

    wf.fc       = cfg.waveform.fc;
    wf.B        = cfg.waveform.B;
    wf.Tc       = cfg.waveform.Tc;
    wf.Fs       = cfg.waveform.Fs;
    wf.N_chirps = cfg.waveform.N_chirps;

    wf.S     = wf.B / wf.Tc;
    wf.Ts    = 1 / wf.Fs;
    wf.N_fast = round(wf.Tc * wf.Fs);
    if abs(wf.Tc * wf.Fs - wf.N_fast) > 1e-6
        error('Tc*Fs = %.4f 非整数, N_fast 采样数不合法', wf.Tc * wf.Fs);
    end
    wf.t_fast = (0:wf.N_fast-1)' * wf.Ts;
    wf.lambda = c / wf.fc;

    % 噪底(式 E8): 各站 NF 应一致, 取第一站
    wf.noise_dbm = -174 + 10*log10(wf.Fs) + cfg.radar(1).NF_db;

    % 帧周期与帧率(帧间空闲由目标帧率推导, Echo_Simulation_Design.md §9.1)
    T_acq = cfg.radar(1).ntx * wf.N_chirps * wf.Tc;        % 采集时长(背靠背帧周期)
    if isfield(cfg.waveform, 'frame_rate') && cfg.waveform.frame_rate > 0
        wf.T_idle = max(0, 1/cfg.waveform.frame_rate - T_acq);
        if wf.T_idle == 0
            warning('目标帧率 %.1f Hz 超过采集能力 %.1f Hz, 退化为背靠背', ...
                    cfg.waveform.frame_rate, 1/T_acq);
        end
    else
        wf.T_idle = 0;
    end
    wf.T_frame    = T_acq + wf.T_idle;
    wf.frame_rate = 1 / wf.T_frame;

    % 校验: 最大拍频 vs Fs/2(双基总路径 = 2*max_range + 阵列口径余量)
    R_max = max([cfg.radar.max_range]);
    if strcmp(cfg.array.type, 'ULA')
        d_arr = (max([cfg.radar.ntx]) * max([cfg.radar.nrx]) - 1) * max([cfg.radar.d_rx]);
    else
        d_arr = max([cfg.radar.ntx]) * max([cfg.radar.d_tx]) + ...
                max([cfg.radar.nrx]) * max([cfg.radar.d_rx]);
    end
    fb_max = wf.S * (2*R_max + d_arr) / c;
    if fb_max >= wf.Fs / 2
        error('最大拍频 %.2f MHz >= Fs/2 = %.2f MHz, 增大 Fs 或减小 max_range', ...
              fb_max/1e6, wf.Fs/2/1e6);
    end
end
