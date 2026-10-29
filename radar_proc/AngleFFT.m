function [az, el, S2db, ax_x, ax_y] = AngleFFT(Rcube, k, ntx, lambda, d_tx, d_rx, nfft_az, nfft_el)
% 距离峰 bin 的 2D 角度谱与峰内插(处理链第三步)
%
% 输入: Rcube (N_fast, S_slots, nrx) 距离谱; k 距离 bin(1-based)
%       ntx; lambda; d_tx; d_rx - 阵列参数(TDM slot->TX 映射 mod(s-1,ntx)+1)
%       nfft_az / nfft_el - 方位/俯仰维补零长度(cfg.proc), >= ntx/nrx
% 输出: az/el (deg) 站系方位/俯仰估计; S2db (nfft_az, nfft_el) fftshift 谱 (dB)
%       ax_x / ax_y 绘图轴(方向余弦 sin(az)cos(el) 与 sin(el))
%
% 阵列排布: TX 竖臂沿站系 z(测俯仰), RX 横臂沿站系 y(测方位)
%
% 原理(Echo_Simulation_Design.md §6b):
%   TDM 重排 slot = tx + ntx*(chirp-1) 与 MATLAB 列主序 reshape 匹配;
%   chirp 维平均(静止目标相位不变);
%   发射程与接收程的距离都随阵元坐标沿视线投影缩短, 相位步进同为负号:
%     omega_tx = -k0*d_tx*sin(el),         (dim1 = TX 竖臂)
%     omega_rx = -k0*d_rx*sin(az)cos(el)   (dim2 = RX 横臂)
%   补零细化谱线: 不改变物理分辨率(孔径决定), 但峰定位精度由谱线间隔
%   ntx/nrx 点提升到 1/nfft 点量级, 且避免 3/4 点短序列抛物线内插失真。
%   补零后保留抛物线内插进一步提精度。

    if nargin < 7, nfft_az = max(64, 2^nextpow2(ntx)); end
    if nargin < 8, nfft_el = max(64, 2^nextpow2(nrx)); end
    k0 = 2*pi / lambda;
    nrx = size(Rcube, 3);
    n_ch = size(Rcube, 2) / ntx;

    % TDM 重排: (slot, rx) -> (tx, chirp, rx)
    A = reshape(squeeze(Rcube(k, :, :)), ntx, [], nrx);
    A = squeeze(mean(A, 2));                         % (ntx, nrx) chirp 平均, 必须 squeeze

    % 加窗 + 补零 2D FFT(dim1 = TX 竖臂 -> 俯仰, 补 nfft_el; dim2 = RX 横臂 -> 方位, 补 nfft_az)
    %A = A .* hamming(ntx) .* hamming(nrx)';   % 汉明窗(用户停用, 矩形窗)

    S2 = fftshift(fft2(A, nfft_el, nfft_az));    % (nfft_el, nfft_az): 行=俯仰, 列=方位
    M = abs(S2);
    S2db = (20*log10(M / max(M(:)))).';          % 转置回 (nfft_az, nfft_el): 行=方位, 列=俯仰(下游 spec2d 约定不变)

    % 2D 峰 + 两维抛物线内插(补零后谱线密, 内插增益小但无害)
    [~, im] = max(M(:));
    [ix, iy] = ind2sub(size(M), im);   % ix: dim1=TX 竖臂(俯仰), iy: dim2=RX 横臂(方位)
    dx = parabolic_delta(M, ix, 1);    % 俯仰维内插
    dy = parabolic_delta(M, iy, 2);    % 方位维内插

    % fftshift 后零频位置 floor(N/2)+1; bin 偏移 -> 数字频率
    bin_el = (ix - 1 - floor(nfft_el/2)) + dx;
    bin_az = (iy - 1 - floor(nfft_az/2)) + dy;
    omega_tx = 2*pi * bin_el / nfft_el;   % TX 竖臂步进 -> sin(el)
    omega_rx = 2*pi * bin_az / nfft_az;   % RX 横臂步进 -> sin(az)cos(el)

    % 数字频率 -> 角度(步进同为负号, 实测验证见 §6b):
    se  = max(min(-omega_tx / (k0 * d_tx), 0.999), -0.999);  % sin(el)
    sac = max(min(-omega_rx / (k0 * d_rx), 0.999), -0.999);  % sin(az)cos(el)
    el = asind(se);
    az = asind(max(min(sac / sqrt(1 - se^2), 0.999), -0.999));

    % 绘图轴: 方向余弦 = -bin * lambda/(d*N), d = lambda/2 时步长 = 2/N
    % 方位轴 ax_x 由 RX 横臂(d_rx)决定, 俯仰轴 ax_y 由 TX 竖臂(d_tx)决定
    ax_x = -((0:nfft_az-1) - floor(nfft_az/2)) * lambda / (d_rx * nfft_az);
    ax_y = -((0:nfft_el-1) - floor(nfft_el/2)) * lambda / (d_tx * nfft_el);
end

% ================================================================
function d = parabolic_delta(M, i, dim)
% 沿指定维的抛物线内插偏移(幅度谱)
    if dim == 1
        a = M(i-1, :);  b = M(i, :);  c = M(i+1, :);
    else
        a = M(:, i-1).';  b = M(:, i).';  c = M(:, i+1).';
    end
    den = a - 2*b + c;
    if den ~= 0
        d = 0.5 * (a - c) / den;
    else
        d = 0;
    end
    d = max(min(d, 1), -1);
end
