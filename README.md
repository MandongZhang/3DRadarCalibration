# 3DRadarCalibration 双雷达自校准仿真

天花板边缘双雷达(同侧 x 两端布置,L 型阵列,3TX×4RX TDM-MIMO 口径)外参自校准的**信号级**仿真程序,含两种运行模式:静态单帧(`main`)与 Z/C/M 轨迹遍历(`main_motion`)。

- **仿真内容**:场景布局(6D 位姿 + L 型阵列)+ 方体目标可见面散射建模 + 每站独立生成中频回波 cube + 单站 FFT 处理链(距离/角度)+ 运动目标轨迹遍历与逐帧检测。
- **设计文档**:[Echo_Simulation_Design.md](Echo_Simulation_Design.md)(式 E1–E19 推导与代码映射)。
- **理论文档**:`Texfiles/2026_Multi-site FMCW Radar Research/` 下
  Research_Ideas.md(四支柱路线)与 Research_Idea-Calibration.md(式 (1)–(30) 推导)。

## 运行

MATLAB 2025b,无工具箱依赖。

```matlab
cd 'D:\BaiduSyncdisk\Labfile\matfile\Multi-site Radar\3DRadarCalibration'
main          % 静态单帧仿真(六项自检 + 三张诊断图)
main_motion   % 运动仿真: Z/C/M 轨迹遍历(逐帧检测 + tracking 图)
```

或命令行 `matlab -batch "main"`、`matlab -batch "main_motion"`。
静态结果存 `results/`(`<date>_echo.mat`,含 cfg/scene/tgt/wf/echo);
运动结果存 `<date>_motion.mat`(traj/pt_hist/truth_hist/detected/P_hist/fuse 等)
与 `<date>_tracking.png`、`<date>_fusion.png`。

---

## 一、设计思路

### 1.1 信号级建模

逐面元双基近场回波相干叠加,SNR、角扩展、强回波等特性全部从几何与雷达方程自然涌现,不人为注入。方体目标的两条特性由此保证:

- **体积大 → 反射强**:可见面元数与单面元 RCS(PO 平板 ∝ A_f²)随尺寸增长
- **体积大 → 占角范围宽**:目标张角内的面元各自回波,角度谱展宽从几何涌现

### 1.2 场景设定

| 设定 | 程序体现 |
|---|---|
| 天花板边缘布置 | `radar(i).p_true` 同侧 x 两端(y = 0.3 m),波束指向地面中心 |
| 两站同高、互不对射 | 无站间互扰建模,每站独立生成回波 |
| 安装存在偏差 | `d_yaw/d_pitch/d_roll` 真值,`p_prior/R_prior` 偏差置零参考位姿 |
| L 型阵列 | TX 沿站系 z 轴(竖臂,测俯仰)、RX 沿站系 y 轴(横臂,测方位),半波长间隔,ntx/nrx 自由配置 |
| 2D TDM-MIMO 虚拟阵 | 虚拟阵元 = t_p + r_q,方位俯仰都可测(元数据 `virt`) |
| 方体目标 | 轴对齐,支持数组,朝向字段预留;静态(main)或沿 Z/C/M 轨迹遍历(main_motion) |
| mirror 散射相位 | 面元初相全零,config 可切 random |

### 1.3 数据流单向

```
config → scene(位姿+阵列) → targets(面元) → waveform → echo(每站 cube)
      → radar_proc(每站距离/角度) → 自检 → 出图 → 存盘
```

所有随机性集中在 `rng(seed)` 之后(散射相位、噪声),结果可复现。雷达为 N 站数组,加第三站只改 `config.m`。

### 1.4 运动仿真(轨迹遍历)

运动入口 `main_motion`:轨迹看成**多个静态点的遍历**——外层循环逐帧把目标
平移到轨迹点(`Move_target`),每帧内部与静态流程完全一致(回波生成 → 处理),
**循环内不画图**,循环后按 `cfg.sim.plot` 开关统一出图。停-走冻结与多普勒
可忽略论证见设计文档 §9(式 E18/E19)。轨迹形状 Z/C/M 由 `cfg.target.traj.type`
选择,检测判定阈值 `cfg.motion.snr_th`。

### 1.5 双站点云融合自校准(SVD)

逐帧检测 (r, az, el) 经 `Det_to_points` 转各站**站系直角坐标**并逐帧累积成
点云,帧序即两站配对关系(同一帧观测同一物理点)。`Fuse_SVD` 用
Kabsch/Umeyama 闭式解(式 E22)求站 2 → 站 1 的旋转 $\hat R$ 与平移
$\hat t$,真值对照为 $R_1'R_2$ 与 $R_1'(p_2-p_1)$(式 E21);配合站 1 先验
即得绝对外参。实测:534–551 配对帧下旋转误差 0.28–0.45°,平移误差
68–84 mm,详见设计文档 §10。

---

## 二、代码结构

| 文件 | 块 | 输入 | 输出 | 公式/说明 |
|---|---|---|---|---|
| `config.m` | — | — | `cfg` | 全部参数,唯一入口 |
| `main.m` | — | — | 打印 + results | 静态流水线 + 自检 + 存盘 |
| `scene/Scene_radarPoses.m` | 一 | `cfg` | `scene` | 6D 位姿 + L 型阵元房间系坐标,式 E10 |
| `scene/Scene_targets.m` | 一 | `cfg` | `tgt` | 方体 6 面面元网格(位置/法向/面积/相位) |
| `scene/Plot_scene.m` | 一 | `cfg`,`scene`,`tgt` | 图 1 | 房间 + 方体面元按照明着色 + 阵列放大插图 |
| `scene/Scene_trajectory.m` | 一(运动) | `cfg` | `traj` | Z/C/M 锚点 → 弧长匀速采样到帧网格(E18) |
| `scene/Move_target.m` | 一(运动) | `tgt` | 平移后 tgt | pos/facet.pos 平移,phase 不动(相位是目标属性) |
| `signal/Waveform_derive.m` | 二 | `cfg` | `wf` | 斜率/快时间轴/噪底,式 E1 E2 校验 |
| `signal/Antenna_gain.m` | 二 | 方向 | 增益 | 单元方向图 cos^g_p,式 E12 |
| `signal/Echo_generate.m` | 二 | `cfg`,`scene`,`tgt`,`wf` | `echo` | 核心:双基近场回波 cube,式 E3–E13 |
| `signal/Check_echo.m` | 三 | 同 main | 自检 | 六项自检,含无噪精确验证路径 |
| `signal/Plot_antenna.m` | 诊断 | `cfg` | 图 2 | 单元方向图:dB 滚降 + fov 标定 + g_p 对比 |
| `radar_proc/RangeFFT.m` | 四 | `echo(i)` | 距离谱 | 快时间 FFT + 距离轴 + 非相干积累 |
| `radar_proc/AngleFFT.m` | 四 | 距离 bin | az/el | TDM 重排 + 2D 角度 FFT,式 E17 |
| `radar_proc/Process_cube.m` | 四 | `scene`,`echo` | `det` | 纯 FFT 主入口:谱最大峰 → (r, az, el) |
| `radar_proc/Truth_project.m` | 四 | `scene`,`tgt` | 真值 | 目标中心站系投影(对照) |
| `radar_proc/Plot_processing.m` | 四 | `det`,`truth` | 图 3 | 每站距离谱 + 2D 角度谱 |
| `radar_proc/Plot_tracking.m` | 四(运动) | `pt_hist` 等 | 跟踪总图 | 3D 轨迹 + SNR 时间线 + r/az/el 对照 + 检测时间线 |
| `radar_proc/Det_to_points.m` | 四(融合) | `pt_hist` 等 | 站系点云 | (r,az,el) → 站系直角坐标,逐帧累积(式 E20) |
| `radar_proc/Fuse_SVD.m` | 四(融合) | `P1`,`P2` | `fuse` | 配对点云 Kabsch/Umeyama SVD 配准(式 E22) |
| `radar_proc/Roomframe_plane.m` | 四(融合) | 站 1 内点云 | `R1`,`p1` | 平面拟合构造房间系(显示用): z=0 取目标平面, 雷达离面高度 d 可观 |
| `radar_proc/Plot_fusion.m` | 四(融合) | `P_hist` 等 | 融合图 | 点云配准叠加 + 逐帧残差 |
| `main_motion.m` | — | 同 main + `traj` | 逐帧检测记录 + results | 静态点遍历循环(内不画图)+ 汇总 + tracking 图 |
| `Echo_Simulation_Design.md` | 文档 | — | — | 设计文档,公式到代码映射 |

---

## 三、变量分析

### 3.1 `cfg`(config.m,唯一参数入口)

| 字段组 | 字段 = 默认值 | 单位 | 含义 |
|---|---|---|---|
| `room` | length=8, width=6, height=2.8 | m | 房间尺寸 |
| `waveform` | fc=77e9, B=1.5e9, Tc=50e-6, Fs=5e6 | Hz,s | 载频/带宽/chirp 时宽/采样率 |
| | N_chirps=64 | — | 每 TX chirp 数 |
| `radar(i)` | p_true=[0.3,0.3,2.8] / [7.7,0.3,2.8] | m | 真实位置(同侧 x 两端) |
| | d_yaw=0/3, d_pitch=0/2, d_roll=0/0 | deg | 安装偏差真值 |
| | ntx=3, nrx=4, d_tx=d_rx=λ/2 | —, m | L 型阵列规模与间距 |
| | fov=60, g_p=3.32 | deg, — | 方向图支撑半锥角与滚降 |
| | Pt_dbm=12.5, G0_dbi=8, NF_db=15 | dBm, dBi, dB | 雷达方程常数 |
| | max_range=10 | m | 波形校验与诊断用(非硬门限) |
| `target` | facet_size=0.05, phase_mode='mirror' | m, — | 面元尺寸与散射相位模式 |
| | targets(1).pos=[3,5,0.1], .size=[0.1,0.1,0.1] | m | 方体中心与三向尺寸 |
| `target.traj` | type='Z', speed=2, n_frames_cap=600 | —, m/s, — | 轨迹形状/匀速/帧数上限(main_motion 用); 高度联动 targets(1).pos(3), 无独立 z 参数 |
| `motion` | snr_th=12 | dB | 检出判定峰SNR阈值(main_motion 用) |
| `sim` | seed=42, save=true, plot=true | — | 随机种子、存盘开关、画图总开关 |

### 3.2 `scene`(Scene_radarPoses 输出)

| 字段 | 维度 | 含义 |
|---|---|---|
| `radar(i).p_true / R_true` | (1,3) / (3,3) | 真实位姿,R = Rz(yaw)Ry(pitch)Rx(roll),boresight = R(:,1) |
| `radar(i).p_prior / R_prior` | (1,3) / (3,3) | 参考位姿(安装偏差置零) |
| `radar(i).tx_pos / rx_pos` | (3,ntx) / (3,nrx) | 阵元房间系坐标(位姿原点 = L 型拐角) |
| `radar(i).boresight` | (3,1) | 波束指向 |

坐标约定:房间系 W 原点墙角、z 向上;站系 S 前 x 右 y 上 z;`C^s = R'(C − p)`。

### 3.3 `tgt`(Scene_targets 输出)

| 字段 | 维度 | 含义 |
|---|---|---|
| `tgt(k).pos / size / R` | (1,3) | 方体中心/尺寸/朝向(当前 eye(3)) |
| `tgt(k).facet.pos / nrm / area` | (M,3)/(M,3)/(M,1) | 面元位置、外法向、面积 |
| `tgt(k).facet.phase` | (M,1) | 散射初相(mirror 全零) |
| `tgt(k).M` | 标量 | 面元数(默认 0.1 m 方体 6 面 24) |

### 3.4 `echo`(Echo_generate 输出,每站一份)

| 字段 | 维度 | 含义 |
|---|---|---|
| `echo(i).data` | (N_fast, ntx·N_chirps, nrx) | **radar cube**,复双精度,噪底归一化为单位复功率 |
| `echo(i).tx_of_slot` | (S,1) | slot → TX 映射 mod(s−1,ntx)+1 |
| `echo(i).virt` | (3, ntx·nrx) | 虚拟阵元房间系坐标,ch = (p−1)·nrx + q |
| `echo(i).t_fast / fc / S / lambda / Fs` | — | 波形透传 |
| `echo(i).noise_dbm` | dBm | 物理噪底参考(绝对量换算用) |
| `echo(i).n_facet` | (1,ntx) | 每 TX 参与面元数(诊断) |

量测模型:逐面元双基近场 `tau = (R_tx + R_rx)/c`,相干叠加见设计文档式 E4。

### 3.5 `wf`(Waveform_derive 输出)

S = B/Tc = 3e13 Hz/s,N_fast = 250,λ = 3.9 mm,噪底 −92 dBm(带宽取 Fs)。
校验:最大双基拍频 < Fs/2;Tc·Fs 为整数。

### 3.6 `traj`(Scene_trajectory 输出)

| 字段 | 维度 | 含义 |
|---|---|---|
| `traj.type` | — | 轨迹形状('Z'/'C'/'M') |
| `traj.t` | (Nf,1) | 帧时刻 (s),t_n = (n−1)·T_frame,T_frame = ntx·N_chirps·Tc |
| `traj.pos` | (Nf,3) | 方体中心房间系坐标,每帧一个静态点 |
| `traj.wp` | (Wp,3) | 形状锚点(绘图参考) |

### 3.7 `fuse`(Fuse_SVD 输出,存于 motion.mat)

| 字段 | 维度 | 含义 |
|---|---|---|
| `fuse.R` | (3,3) | 站 2 → 站 1 旋转估计(真值 $R_1'R_2$) |
| `fuse.t` | (3,1) | 站 2 → 站 1 平移估计(真值 $R_1'(p_2-p_1)$) |
| `fuse.res` | (N,1) | 逐配对帧残差范数 (m) |
| `fuse.rmse` | — | 配对残差 RMSE (m) |

motion.mat 另存 `P_hist/T_hist (Nf,n_st,3)` 站系点云与真值云、`both` 配对掩码、
`R_rel/t_rel` 相对外参真值。

---

## 四、自检项(Check_echo,main 运行时打印)

| 自检 | 策略 | 容差 |
|---|---|---|
| 1 拍频 bin | 无噪,单面元上 bin,抛物线插值 | 0.05 bin |
| 2 峰相位(含 −Sτ²/2) | 无噪,解析式核对 | 5e-3 rad(实测 ~1e-14) |
| 3 通道间相位差 | 无噪,近场解析;平面波近似作 INFO 对照 | 1e-3 rad(实测 ~1e-13) |
| 4 峰 SNR = E7 + 10log10(N_fast) | 有噪,大面元提 SNR | ±2 dB |
| 5 大目标谱功率增强 | 有噪,2m/1m 对比 vs 解析预测 | ±2 dB |
| 6 谱功率定标 | 非相干预测为基准,近场相干干涉损失计入容差 | ±12 dB |

静态模型在运动仿真的每帧内严格成立,六项自检同样覆盖 main_motion。

## 五、参考

- Research_Ideas.md:支柱 A 距离主导 / B 双 Doppler 采集 / C 离轴-距离误差建模 / D 场景适配
- Research_Idea-Calibration.md:式 (1)–(30) 推导与流程
- Echo_Simulation_Design.md:信号级模型 E1–E19 与代码映射
- PaperReading/2022_EuRAD_Grebner:σ 框架与嵌套优化来源
