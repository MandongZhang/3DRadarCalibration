function out = Register_iterSVD(D1, D2, R1p, p1, R2p, p2, h, max_iter)
% ULA 平面约束迭代配准(新算法核心)
%
% 模型: 两站 ULA 观测 (r, u_y); 目标轨迹在房间系水平面 z = h 上(物理前提);
%       平面在站 i 系的表达 (nu_i, d_i) 依赖站姿态, 与相对位姿耦合,
%       用定点迭代解: 圆弧求交 → Fuse_SVD → 平面重拟合, 初值取安装先验。
%       站 1 为参考站: 其姿态先验误差不参与修正(与 ref 的相对位姿口径一致)。
%
% 输入: D1/D2 (N,2) [r, u_y] 双站同时检出帧(帧序即配对);
%       R1p/p1, R2p/p2 两站先验位姿(位置无先验误差, 仅姿态有偏差);
%       h 轨迹平面高度(房间系); max_iter 最大迭代轮数(默认 8)
% 输出: out.R/out.t 相对位姿估计(站2->站1, 同 Fuse_SVD 约定);
%       out.res/out.rmse 末轮残差; out.n_iter/out.conv;
%       out.nu1/out.d1 平面(站1系); out.Cs1/out.Cs2 末轮交点云(站系);
%       out.n_inlier 3σ 内点数
    if nargin < 8, max_iter = 8; end
    R1 = R1p;  R2 = R2p;
    % 平面初值(站1系): 房间 z=h 平面在先验姿态下的表达
    nu1 = R1p(3,:);  d1 = h - p1(3);
    R_old = eye(3);  t_old = zeros(3,1);
    f = [];  Cs1 = [];  Cs2 = [];
    conv = false;
    for it = 1:max_iter
        % 当前相对位姿估计下的平面(站2系)
        R_rel = R1.'*R2;   t_rel = R1.'*(p2 - p1).';
        nu2 = nu1*R_rel;   d2 = d1 - t_rel.'*nu1.';
        % 圆弧求交(各站自身系)
        Cs1 = RayPlanePoints(nu1, d1, D1);
        Cs2 = RayPlanePoints(nu2, d2, D2);
        % 刚体配准(站2 -> 站1, ref 的 Fuse_SVD 原样复用)
        f = Fuse_SVD(Cs1, Cs2);
        dR = norm(f.R - R_old, 'fro');   dt = norm(f.t - t_old);
        R_old = f.R;  t_old = f.t;
        conv = (dR < 1e-9 && dt < 1e-9);
        % 平面重拟合: 站2 云经新位姿映射到站1系, 与站1云联合 TLS 拟合
        C2in1 = Cs2*f.R.' + f.t.';
        [nu1, d1] = tls_plane([Cs1; C2in1]);
        R2 = R1*f.R;
        if conv, break; end
    end
    out.R = f.R;  out.t = f.t;  out.res = f.res;  out.rmse = f.rmse;
    out.n_iter = it;  out.conv = conv;
    out.plane1 = struct('nu', nu1, 'd', d1);
    out.Cs1 = Cs1;  out.Cs2 = Cs2;
    % 3σ 剔除后重解一轮(与 ref 的两遍法对应)
    inl = out.res < 3*out.rmse;
    out.n_inlier = nnz(inl);
    if nnz(inl) > 10 && nnz(inl) < numel(inl)
        nu2 = out.plane1.nu*out.R;
        d2 = out.plane1.d - out.t.'*out.plane1.nu.';
        Cs1 = RayPlanePoints(out.plane1.nu, out.plane1.d, D1(inl,:));
        Cs2 = RayPlanePoints(nu2, d2, D2(inl,:));
        f = Fuse_SVD(Cs1, Cs2);
        out.R = f.R;  out.t = f.t;  out.res = f.res;  out.rmse = f.rmse;
        out.Cs1 = Cs1;  out.Cs2 = Cs2;
    end
end

function [nu, d] = tls_plane(P)
% TLS 平面拟合: 过质心, 最小奇异向量为法向(行向量)
    mu = mean(P, 1);
    X = P - mu;
    [~, ~, V] = svd(X, 0);
    nu = V(:, end).';
    d = mu*nu.';
end
