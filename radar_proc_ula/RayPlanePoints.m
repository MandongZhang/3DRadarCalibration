function [Cs, ok] = RayPlanePoints(nu, d, Dm)
% ULA 测量圆弧 × 平面求交(平面在测量的站系内表达: C_s * nu' = d)
%
% 输入: nu (1,3) 平面单位法向(站系); d 平面偏距; Dm (N,2) [r, u_y]
% 输出: Cs (N,3) 站系三维交点( = r*u ); ok (N,1) 交点存在且取到分支
%
% 原理: 站系方向 u = [c*cosγ; u_y; c*sinγ], c = sqrt(1-u_y^2);
%       平面约束 r*(u·nu') = d → nu(1)c cosγ + nu(3)c sinγ = d/r - nu(2)u_y
%       两支解取前向分量 u(1) 较大者(目标位于 boresight ±90° 内)
    N = size(Dm, 1);
    Cs = nan(N,3);  ok = false(N,1);
    r = Dm(:,1);  uy = Dm(:,2);
    c = sqrt(max(1 - uy.^2, 0));
    A = nu(1)*c;
    B = nu(3)*c;
    Ct = d./r - nu(2)*uy;
    AB = sqrt(A.^2 + B.^2);
    has = (AB > 0) & (abs(Ct) <= AB);
    ac = acos(max(min(Ct./max(AB, eps), 1), -1));
    g0 = atan2(B, A);
    for n = 1:N
        if ~has(n), continue; end
        g1 = g0(n) + ac(n);  g2 = g0(n) - ac(n);
        u1 = [c(n)*cos(g1); uy(n); c(n)*sin(g1)];
        u2 = [c(n)*cos(g2); uy(n); c(n)*sin(g2)];
        if u1(1) >= u2(1), uu = u1; else, uu = u2; end
        Cs(n,:) = (r(n)*uu).';
        ok(n) = true;
    end
end
