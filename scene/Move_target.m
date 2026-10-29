function tgt = Move_target(tgt, dpos)
% 平移目标(轨迹遍历用): pos 与 facet.pos 整体平移 dpos, 其余属性不动
%
% 等效于修改 cfg.target.targets(k).pos 后重跑 Scene_targets, 但保留
% "散射相位是目标属性"的约定(random 相位模式不逐帧重抽), 且免去面元重建。
% 输入: tgt - Scene_targets 输出; dpos (1,3) 平移量 (房间系, m)
% 输出: tgt - 平移后的目标结构体数组

    for k = 1:numel(tgt)
        tgt(k).pos = tgt(k).pos + dpos;
        tgt(k).facet.pos = tgt(k).facet.pos + repmat(dpos, size(tgt(k).facet.pos, 1), 1);
    end
end
