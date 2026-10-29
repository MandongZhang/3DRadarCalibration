function Compare_ref_ula(date)
% ref(L 型闭式 SVD) 与 ULA(平面约束迭代) 配准精度并排对照
% 用法: Compare_ref_ula('2026-10-29') 或 Compare_ref_ula()(自动取两类结果共有的最新日期)
    if nargin < 1
        dr = dir(fullfile('results', '*_motion.mat'));
        du = dir(fullfile('results', '*_motion_ula.mat'));
        dates = intersect({dr.name}, {du.name});
        assert(~isempty(dates), 'results/ 中没有同日期的 ref 与 ULA 结果');
        date = dates{end}(1:10);
    end
    A = load(fullfile('results', [date '_motion.mat']));
    B = load(fullfile('results', [date '_motion_ula.mat']));
    S = {A, B};
    names = {'ref (L闭式SVD)', 'ULA (平面约束迭代)'};
    fprintf('==== 配准精度对照 (%s) ====\n', date);
    fprintf('%-22s %14s %14s %16s %10s\n', '算法', '旋转误差(deg)', '平移误差(mm)', '内点RMSE(mm)', '配对帧');
    for k = 1:2
        Sk = S{k};
        eR = acosd(max(min((trace(Sk.R_rel.' * Sk.fuse.R) - 1)/2, 1), -1));
        eT = norm(Sk.fuse.t - Sk.t_rel) * 1e3;
        rmse = Sk.fuse.rmse * 1e3;
        np = numel(Sk.fuse.res);
        fprintf('%-22s %14.4f %14.2f %16.2f %10d\n', names{k}, eR, eT, rmse, np);
    end
end
