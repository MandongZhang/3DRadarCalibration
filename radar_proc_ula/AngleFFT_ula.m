function [uy, S1db, ax] = AngleFFT_ula(Rcube, k, ntx, lambda, d_v, nfft)
% ULA 1D 角度处理: 虚拟阵(沿站系 y, 间距 d_v = λ/2)1D FFT → u_y = sin(az)cos(el)
%
% 原始观测量是方向余弦 u_y, 不做 az = asind(u_y/cos(el)) 反解
% (el 未知, 直接反解会引入 ~2-4° 系统偏差; el 由配准迭代的平面求交找回)
%
% 输入: Rcube (N_fast, ntx*N_chirps, nrx) 距离谱; k 距离 bin(1-based)
%       ntx; lambda; d_v 虚拟阵元间距(= d_rx); nfft 补零长度
% 输出: uy 峰方向余弦(钳制 ±1); S1db (nfft,1) fftshift 谱(dB); ax (nfft,1) u 轴
%
% 相位步进符号与 E17 一致: 路程随阵元坐标沿视线缩短, omega = -k0*d_v*u_y
    if nargin < 6, nfft = 256; end
    k0 = 2*pi/lambda;
    nrx = size(Rcube, 3);
    % TDM 重排: (slot, rx) -> (tx, chirp, rx); chirp 平均(必须 squeeze)
    A = reshape(squeeze(Rcube(k,:,:)), ntx, [], nrx);
    A = squeeze(mean(A, 2));                       % (ntx, nrx)
    % 虚拟阵排序: ch = (p-1)*nrx + q <-> 虚拟位置 ((p-1)*nrx + q-1)*d_v
    A1d = reshape(A.', 1, []);                     % (1, ntx*nrx), 下标即通道号
    S1 = fftshift(fft(A1d.', nfft));
    M = abs(S1);
    S1db = 20*log10(M/max(M));
    [~, im] = max(M);
    if im > 1 && im < nfft
        den = M(im-1) - 2*M(im) + M(im+1);
        if den ~= 0, dl = 0.5*(M(im-1)-M(im+1))/den; else, dl = 0; end
        dl = max(min(dl,1),-1);
    else
        dl = 0;
    end
    bin = (im - 1 - floor(nfft/2)) + dl;
    uy = max(min(-2*pi*bin/nfft/(k0*d_v), 1), -1);
    ax = -((0:nfft-1) - floor(nfft/2)) * lambda/(d_v*nfft);
end
