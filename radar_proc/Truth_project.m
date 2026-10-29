function T = Truth_project(scene, tgt)
% 目标中心真值投影到各站视角(处理结果对照用)
%
% 输入: scene - Scene_radarPoses; tgt - Scene_targets
% 输出: T(i).truth (K,3) [r (m) | az (deg) | el (deg)], 站 i 视角下各目标中心
%
% 投影(站系): C^s = R'(C - p)', r = ||C^s||,
% az = atan2(y^s, x^s), el = atan2(z^s, hypot(x^s, y^s))

    K = numel(tgt);
    for i = 1:numel(scene.radar)
        rad = scene.radar(i);
        T(i).truth = nan(K, 3);
        for k = 1:K
            Cs = rad.R_true' * (tgt(k).pos.' - rad.p_true.');
            r  = norm(Cs);
            az = atan2d(Cs(2), Cs(1));
            el = atan2d(Cs(3), hypot(Cs(1), Cs(2)));
            T(i).truth(k, :) = [r, az, el];
        end
    end
end
