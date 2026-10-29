function [Rcube, r_ax, P_acc] = RangeFFT(echo_i)
% 快时间维 FFT + 距离轴 + 非相干积累谱(处理链第一步)
%
% 输入: echo_i - echo(i), Echo_generate 输出
% 输出: Rcube (N_fast, S_slots, nrx) 距离谱(复)
%       r_ax  (N_fast,1) 距离轴 (m): r_k = k*Fs/N_fast * c/(2S)
%       P_acc (N_fast,1) 非相干积累谱(线性功率): mean over (slot, rx) |X|^2
%
% 距离换算用单站近似 r = c*fb/(2S); 双基差(站内收发间距毫米级)可忽略。

    c = 299792458;
    X = fft(echo_i.data, [], 1);
    N_fast = size(X, 1);
    r_ax = (0:N_fast-1).' * (echo_i.Fs / N_fast) * c / (2 * echo_i.S);
    P_acc = mean(abs(X).^2, [2 3]);
    Rcube = X;
end
