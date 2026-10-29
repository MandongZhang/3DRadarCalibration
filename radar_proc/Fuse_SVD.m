function out = Fuse_SVD(P1, P2)
% 配对点云 SVD 刚体配准(Kabsch/Umeyama, 无尺度): 求 R, t 使 R*P2 + t ≈ P1
%
% 输入: P1, P2 (N,3) 逐帧配对的两站站系点(同一物理点在站 1 / 站 2 系的坐标,
%       帧序即对应关系, 免 ICP 最近邻迭代)。无效帧(未检出)须在调用前剔除。
% 输出: out.R    (3,3) 旋转, 站 2 系 -> 站 1 系
%       out.t    (3,1) 平移, 站 2 系原点在站 1 系的坐标
%       out.res  (N,1) 各配对残差范数 (m)
%       out.rmse 配对残差 RMSE (m)
%
% 解法(Echo_Simulation_Design.md §10): 去质心 X = P2 - mean(P2), Y = P1 - mean(P1)
%   H = X'Y;  SVD H = U S V';
%   R = V diag(1,1,sign det(VU')) U';   t = mean(P1)' - R * mean(P2)'

    mu1 = mean(P1, 1).';                   % (3,1) 站 1 云质心
    mu2 = mean(P2, 1).';                   % (3,1) 站 2 云质心
    H   = (P2 - mu2.')' * (P1 - mu1.');    % (3,3) 互协方差 H = X'Y(Kabsch)
    [U, ~, V] = svd(H, 0);
    D  = diag([1, 1, sign(det(V * U.'))]); % 反射校正
    out.R   = V * D * U.';
    out.t   = mu1 - out.R * mu2;
    out.res = vecnorm(P1 - (P2 * out.R.' + out.t.'), 2, 2);
    out.rmse = sqrt(mean(out.res.^2));
end
