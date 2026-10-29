# 信号级回波仿真设计(公式推导 → 代码映射)

本文是信号级仿真的设计文档:每台雷达根据对目标的照射情况独立生成中频回波
(radar cube),处理链(range/angle FFT)在 `radar_proc/` 中实现,运动目标的
轨迹遍历见 §9。

场景:8×6×2.8 m 房间,天花板边缘双雷达(同侧 x 两端),L 型阵列(TDM-MIMO),
方体目标(可见面散射中心),安装偏差真值 + 偏差置零参考位姿,随机种子可复现。

## 0. 记号约定与符号对照

### 0.1 下标与上标

| 记号 | 含义 | 取值 |
|---|---|---|
| $i$ | 雷达站编号 | $1..N$,本项目 $N=2$ |
| $p,\ q$ | TX / RX 阵元编号 | $1..n_{tx}$,$1..n_{rx}$ |
| $m$ | 面元编号 | $1..M$ |
| $k$ | 目标(方体)编号 | $1..K$,默认 1 |
| $s$ | TDM slot 编号 | $1..n_{tx}N_{chirps}$ |
| $W,\ S$ | 房间系 / 站系 | 站系前 $x$(波束)右 $y$ 上 $z$ |
| $u$ | 单位方向向量 | 视线方向,归一化 |

### 0.2 向量与矩阵实体

| 符号 | 维度 | 坐标系 | 含义 |
|---|---|---|---|
| $p_i$ | 3×1 | 房间系 | 站 $i$ 位置(阵列相位中心 = 阵中心) |
| $R_i$ | 3×3 | 站系→房间系 | 旋转,$R = R_z(\text{yaw})R_y(\text{pitch})R_x(\text{roll})$ |
| $p^{tx}_p$ | 3×1 | 房间系 | TX 阵元 $p$ 位置 |
| $p^{rx}_q$ | 3×1 | 房间系 | RX 阵元 $q$ 位置 |
| $C_k$ | 3×1 | 房间系 | 方体 $k$ 中心 |
| $p_m$ | 3×1 | 房间系 | 面元 $m$ 中心 |
| $n_m$ | 3×1 | 房间系 | 面元外法向(单位) |
| $A_f$ | 标量 | — | 面元面积 (m²) |
| $\psi_m$ | 标量 | — | 面元散射初相(mirror 模式全 0) |
| $\sigma_m$ | 标量 | m² | 面元 RCS |
| $u_{tx},\ u_{rx}$ | 3×1 | — | 面元到 TX / RX 的单位方向 |
| $\tau_{mpq}$ | 标量 | s | 双基双程延迟 |
| $s_{IF}(t)$ | 复 | — | 中频信号,快时间 $t$ |

### 0.3 公式符号与代码对照

| 公式符号 | 代码变量 | 产出函数 |
|---|---|---|
| $p^{tx}_p,\ p^{rx}_q$ | `scene.radar(i).tx_pos / rx_pos` (3,ntx)/(3,nrx) | Scene_radarPoses |
| $p_m,\ n_m,\ A_f,\ \psi_m$ | `tgt(k).facet.pos/nrm/area/phase` | Scene_targets |
| $S,\ t,\ \lambda$ | `wf.S / t_fast / lambda` | Waveform_derive |
| $G(u)$ | `Antenna_gain` 输出 | Antenna_gain |
| $s_{IF}$ 采样阵 | `echo(i).data` (N_fast, ntx·N_chirps, nrx) | Echo_generate |
| slot→TX 映射 | `echo(i).tx_of_slot` | Echo_generate |
| 虚拟阵元 $v_{pq}$ | `echo(i).virt` (3, ntx·nrx) | Echo_generate |
| 物理噪底 | `wf.noise_dbm` / `echo(i).noise_dbm` | Waveform_derive |

## 1. FMCW 波形模型

TX chirp(第 $p$ 个 TX,快时间 $t \in [0, T_c)$):

$$
s_T(t) = \exp\Big(j2\pi\big(f_c\,t + \tfrac{S}{2}t^2\big)\Big),\qquad
S = \frac{B}{T_c} \tag{E1}
$$

参数默认值(fc = 77 GHz,B = 1.5 GHz,Tc = 50 µs,Fs = 5 MHz,N_chirps = 64):

$$
S = 3\times10^{13}\ \text{Hz/s},\qquad
N_{fast} = T_c F_s = 250,\qquad
\Delta R = \frac{c}{2B} = 10\ \text{cm}
$$

快时间轴 $t = (0:N_{fast}-1)\,T_s$,$T_s = 1/F_s$。

波形合法性校验(失败即 error):最大拍频须落在采样率内,

$$
f_{b,\max} = S\cdot\frac{2\,R_{max} + D_{arr}}{c} < \frac{F_s}{2} \tag{E2}
$$

其中 $D_{arr}$ 为阵列口径余量。双基总路径取代单程距离参与校验。

## 2. 双基几何与中频信号

### 2.1 双基延迟(近场精确)

面元 $m$、TX $p$、RX $q$ 的双程路径分段计算:

$$
R_{tx} = \|p_m - p^{tx}_p\|,\qquad
R_{rx} = \|p^{rx}_q - p_m\|,\qquad
\tau_{mpq} = \frac{R_{tx} + R_{rx}}{c} \tag{E3}
$$

站内 TX/RX 相距毫米级,面元距离米级,按近场精确计算,不引入平面波近似。
远场极限下 $\tau$ 退化为 $\big(2R + (t_p + r_q)\cdot u\big)/c$ 的大标形式,
虚拟阵元概念(§5)由此成立。

### 2.2 中频信号

RX 端与本振混频后,残余项为常数相位与拍频线性项:

$$
s_{IF}(t) = \sum_{m \in \mathcal{V}} a_{mpq}\,
\exp\Big(j2\pi\Big(f_c\tau_m + S\tau_m t - \frac{S}{2}\tau_m^2\Big)\Big) \tag{E4}
$$

- 拍频 $f_b = S\tau_m$
- $-S\tau_m^2/2$ 为残余相位项:4 m 路径约 0.045 rad,自检 2 对该项敏感
- $\mathcal{V}$ 为通过可见性与视场判据(§4.3)的面元集合
- 复幅度 $a_{mpq}$ 由 §3.2 雷达方程给出

### 2.3 静止目标的慢时间结构

每帧内目标冻结(运动目标的逐帧处理见 §9),同一 TX 的各 chirp 确定部分完全
相同。故信号部分对每 TX 只合成一次
$(N_{fast} \times n_{rx})$,复制到 chirp 维;噪声逐 slot 独立生成:

$$
\text{data}(:,\,\text{cols},\,:) = \text{repmat}(s_0,\ 1,\ N_{chirps},\ 1) + n(t) \tag{E5}
$$

## 3. 面元散射模型

### 3.1 方体离散化与物理光学 RCS

方体 $k$ 六面网格化,每边 $n = \mathrm{round}(L/\Delta)$ 等分,
$d_{eff} = L/n$ 精确铺满。每个面元视作小平板散射体,物理光学(PO)近似:

$$
\sigma_m = \frac{4\pi}{\lambda^2}A_f^2\cos^2\theta_b,\qquad
\cos\theta_b = \frac{n_m\cdot u_{tx} + n_m\cdot u_{rx}}{\|u_{tx} + u_{rx}\|},\ \ \ge 0 \tag{E6}
$$

$\theta_b$ 为面元法向与 TX/RX 平分线方向的夹角(双基站内近似,站内 TX/RX
视角差 < 0.3°)。mirror 模式散射系数为正实数,面元间干涉相位完全来自传播相位。

**收敛性**:镜面峰 $\sigma_{peak} = 4\pi A_{total}^2/\lambda^2$ 只依赖总面积,
与面元网格无关(自检 6 验证)。**栅格化已知特性**:偏离镜面方向存在角宽约
$\lambda/d$、比峰低约 20 dB 的本底,"大目标占角范围更宽"由此自然涌现,
是模型特性而非缺陷。

config 预留 `phase_mode='random'`:面元初相按 $U[0,2\pi)$ 种子生成,模拟粗糙
目标,总 RCS 趋近非相干叠加。**相位是目标属性,在 Scene_targets 按种子一次
生成,不得逐站重新生成**。

### 3.2 双基雷达方程与噪底定标

$$
P_{rx} = \frac{P_t\,G_{tx}(u)\,G_{rx}(u)\,\lambda^2\,\sigma_m}{(4\pi)^3\,R_{tx}^2\,R_{rx}^2} \tag{E7}
$$

$G$ 逐阵元方向图(§4.2),不是逐站。噪底定标链四步:

1. 物理噪底(kT₀ = −174 dBm/Hz @ 290 K,带宽取 $F_s$,$N_F$ 噪声系数):

$$
P_{noise,\mathrm{dBm}} = -174 + 10\log_{10}(F_s) + N_F \tag{E8}
$$

2. 雷达方程 (E7) 逐 $(m,p,q)$ 算 $P_{rx}$
3. 数据归一化:噪声复功率为 1,$n = (\mathcal{N}_1 + j\mathcal{N}_2)/\sqrt{2}$;
   信号复幅度

$$
a_{mpq} = \sqrt{\frac{P_{rx}}{P_{noise}}}\;e^{j\psi_m},\qquad \text{逐样本 SNR} = |a|^2 \tag{E9}
$$

4. echo 元数据存 `noise_dbm`,下游换算绝对量

**数量级预演**:77 GHz 下 1 m² 平板法向 $\sigma = 4\pi A^2/\lambda^2 \approx$
59 dBsm,$P_t$ = 12.5 dBm、$G$ = 8 dBi ×2、$R$ = 4 m 时 $P_{rx} \approx$ −18 dBm,
噪底约 −92 dBm。镜面方向 SNR 高是物理正确的;非镜面面元由 $\cos^2\theta_b$
压回动态范围。double 精度无饱和问题。

### 4. 天线方向图、可见性与视场

### 4.1 位姿旋转与阵列几何(L 型,TDM-MIMO)

#### 4.1.1 旋转约定与欧拉角(euler2R)

旋转矩阵采用 ZYX 内旋顺序的 yaw-pitch-roll 航空安装约定,矩阵作用自右向
左(向量先滚转、再俯仰、再偏航):

$$
R = R_z(\text{yaw})\,R_y(\text{pitch})\,R_x(\text{roll}) \tag{E14}
$$

$$
R_z = \begin{bmatrix} c_y & -s_y & 0 \\ s_y & c_y & 0 \\ 0 & 0 & 1 \end{bmatrix},\quad
R_y = \begin{bmatrix} c_p & 0 & s_p \\ 0 & 1 & 0 \\ -s_p & 0 & c_p \end{bmatrix},\quad
R_x = \begin{bmatrix} 1 & 0 & 0 \\ 0 & c_r & -s_r \\ 0 & s_r & c_r \end{bmatrix}
$$

其中 $c_y = \cos(\text{yaw})$ 等均为角度制(代码 `cosd/sind`)。选 ZYX 顺序的
三条理由:匹配物理安装顺序(先定水平朝向、再调下俯、最后微调横滚,偏差
真值 $d_{yaw}, d_{pitch}, d_{roll}$ 直接加在对应角上);pitch 与 yaw 解耦(见下);
奇异点(万向锁)在 pitch = ±90°,本场景俯仰约 30° 远离。

**波束指向与 pitch 解耦**(第一列展开):

$$
R(:,1) = R_z R_y \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
= \begin{bmatrix} c_p c_y \\ c_p s_y \\ -s_p \end{bmatrix} \tag{E15}
$$

第三分量只含 pitch,$-s_p$ 给出约定:**pitch > 0 为下视**。yaw 与 pitch 解耦,
故"理想朝向 = 指向地面中心"可两步独立反解(Scene_radarPoses 中
$\boldsymbol{d} = [L/2 - p_x,\ W/2 - p_y,\ -p_z]^{\mathsf T}$):

$$
\text{ideal\_yaw} = \operatorname{atan2}(d_y,\ d_x),\qquad
\text{ideal\_pitch} = -\arcsin\!\Big(\frac{d_z}{\|\boldsymbol{d}\|}\Big) \tag{E16}
$$

真实位姿 = 理想 + 安装偏差:$\text{yaw}_t = \text{ideal\_yaw} + d_{yaw}$ 等;
参考位姿 = 偏差置零。

**roll 的几何地位**:roll 不进第一列(式 E15 与 $c_r, s_r$ 无关)。绕波束轴
滚转不改变波束指向,只转动 L 型竖臂(RX 排布方向),故互不对射场景对绕
基线滚转的敏感度低。

#### 4.1.2 阵列几何

位姿原点取 L 型拐角(TX1 位置)。TX 沿站系 $z$ 轴正向(竖臂,测俯仰)从原点
排布、RX 沿站系 $y$ 轴正向(横臂,测方位)偏移一个阵元起排(拐角处 TX1 与
RX1 间隔 $d_{rx}$,收发不重叠):

$$
t_p = (p-1)\,d_{tx},\qquad
r_q = q\,d_{rx},\qquad
d_{tx} = d_{rx} = \frac{\lambda}{2} \tag{E10}
$$

两臂整体平移不改变阵元相对位置,虚拟阵测角能力不变。TDM 时隙:
slot $s$ 的发射 TX $= \bmod(s-1,\ n_{tx}) + 1$,slot 总数 $n_{tx}N_{chirps}$。
虚拟阵元(元数据,$v_{pq} = t_p + r_q$,通道排序 $\text{ch} = (p-1)n_{rx} + q$):

$$
v_{pq} = t_p + r_q \tag{E11}
$$

远场下 TX 维控俯仰、RX 维控方位,2D 虚拟阵两维都可测。

### 4.2 单元方向图

$$
G(u) = \cos^{g_p}(\psi),\qquad \psi = \arccos\big(u\cdot \hat{b}\big),\quad u = \frac{p_m - p^{ant}}{\|p_m - p^{ant}\|} \tag{E12}
$$

$\hat{b} = R(:,1)$ 为 boresight。$g_p = 3.32$ 沿用(标定到 −10 dB @ 60°)。
fov 为方向图支撑半锥角(默认 60°),超出置 0(软限)。**方向约定**:$u$ 取
阵元指向面元方向;可见性判据 E13 的 $u_{tx},u_{rx}$ 取面元指向阵元方向,
两者互为反向量,勿混用。

### 4.3 可见性判据

面元 $m$ 对 $(p,q)$ 参与叠加须同时满足:

$$
n_m\cdot u_{tx} > 0 \ \text{且}\ n_m\cdot u_{rx} > 0 \ \text{且}\ G\ne 0,\qquad u_{tx} = \frac{p^{tx}_p - p_m}{\|p^{tx}_p - p_m\|},\ \ u_{rx} = \frac{p^{rx}_q - p_m}{\|p^{rx}_q - p_m\|} \tag{E13}
$$

$u_{tx},u_{rx}$ 均由**面元指向阵元**,法向朝阵元即为被照明/可见。掠射面元
不硬删,$\cos^2\theta_b$ 自然压零,硬判据只为省算力。

## 5. Echo_generate 代码映射与向量化

### 5.1 输出结构

```
echo(i).data       (N_fast, ntx*N_chirps, nrx) 复双精度 IF 采样, 噪底归一化
echo(i).tx_of_slot (ntx*N_chirps,1)  slot->TX 映射
echo(i).virt       (3, ntx*nrx)      虚拟阵元房间系坐标
echo(i).t_fast/.fc/.S/.lambda/.Fs    波形透传
echo(i).noise_dbm  物理噪底
echo(i).n_facet    (1,ntx)           每 TX 参与叠加的面元数(诊断)
```

### 5.2 向量化伪代码(循环只到 站 × TX 两层)

```matlab
P, Nrm, Af, phase_f = 全目标面元聚合            % (3,M) (3,M) (M,1) (M,1)
for i = 1:站数
  data = complex(zeros(N_fast, ntx*N_chirps, nrx))
  for p = 1:ntx
    dTx = tx_pos(:,p) − P;  Rtx = vecnorm;  uTx = dTx./Rtx        % (3,M) (M,1) 面元->TX
    dRx = rx_pos − P 扩展 (3,M,nrx);  Rrx, uRx 同理               % (M,nrx) 面元->RX
    vis  = (n·uTx > 0) & (n·uRx > 0)                              % (M,nrx) E13
    Gtx, Grx = Antenna_gain(−uTx, −uRx)                           % 方向图取阵元->面元
    cosb = 平分线夹角(E6), clamp ≥ 0
    sigma = (4π/λ²)·Af²·cosb²                                     % (M,nrx)
    Prx   = Pt·Gtx·Grx·λ²·sigma ./ ((4π)³·Rtx²·Rrx²)              % (M,nrx)
    W     = sqrt(Prx) · vis · exp(j·phase_f)                      % 复幅度 (M,nrx)
    TAU   = (Rtx + Rrx)/c                                         % (M,nrx) 近场精确
    beta  = exp(j2π(fc·TAU − S/2·TAU²))                           % 常相位残余项
    fb    = S·TAU
    s0 = zeros(N_fast, nrx)
    for 块 of 512 面元:                                            % 临时量 ~8 MB
      Eb = exp(j2π·fb(块)·t_fast)   广播 (块, N_fast, nrx)
      s0 += sum(Eb .* (W·beta)(块), 1)
    end
    cols = p : ntx : ntx*N_chirps                                 % TDM 槽位
    noise = (randn + 1i·randn)/√2, (N_fast, N_chirps, nrx)        % 逐 slot 独立
    data(:,cols,:) = repmat(reshape(s0,N,1,nrx), 1, N_chirps, 1) + noise  % E5 先升维再复制
  end
  echo(i).data = data  附元数据
end
```

要点:
- `beta`、`W`、`sigma`、`vis` 均为 M×nrx 小矩阵,廉价
- 3D 临时量只有分块 `Eb`,Mblk 控制峰值内存
- chirp 维零循环(静止目标)
- `fc·TAU` 相位幅值 ~3e3 rad,double 误差 ~1e-12 rad,无需 mod 包装

### 5.3 内存与耗时

| 对象 | 维度 | 体积 |
|---|---|---|
| echo(i).data(默认参数) | 250×192×4 复 double | 3.1 MB/站,两站 6.2 MB |
| 分块 Eb 临时量 | 512×250×4 复 double | 8.2 MB 峰值 |
| results/*.mat | 全部结构 + 元数据 | < 20 MB |

复指数约 4e6 次每站,总耗时秒级。极端放大(N_fast=1024, 128 chirps, 8×8)
→ 134 MB/站,仍可行。

## 6. 自检设计(Check_echo)

无处理链前提下闭环验证回波正确性。自检内部构造覆盖 cfg(单面元/双阵元/上 bin
距离)独立调 Echo_generate,不污染正式数据。上 bin 消除频谱泄漏。

**策略**:自检 1/2/3 用无噪回波(Echo_generate 第 5 参 false)精确验证确定性
模型,实测误差在数值精度(~1e-13 rad);自检 4/5/6 用有噪回波验证统计定标,
容差按噪声涨落放大。

| # | 内容 | 要点 | 容差 |
|---|---|---|---|
| 1 | 快时间 FFT 峰位置 = 理论拍频 bin | 单点单站,距离上 bin,抛物线插值 | < 0.05 bin |
| 2 | 峰相位 = $2\pi(f_c\tau - S\tau^2/2)$ | 无噪上 bin 无 Dirichlet 项;$\tau = 2D/c$ | ±5e-3 rad |
| 3 | 相邻 RX / TX 通道相位差 | 无噪;近场解析 $2\pi(f_c\Delta\tau - S\Delta(\tau^2)/2)$ 主判据;平面波近似(发射程负号、接收程正号)作 INFO 对照 | ±1e-3 rad |
| 4 | 峰 SNR = 雷达方程 + $10\log_{10}N_{fast}$ | 有噪,大面元(1e-3 m²)提 SNR 压涨落 | ±2 dB |
| 5 | 大目标谱功率增强 | 有噪,2m/1m 方体(正面同距)对比 vs 解析预测 | ±2 dB |
| 6 | 谱功率定标 | 非相干逐面元求和为基准;近场方体相邻面元拍频差 < bin 宽时干涉保留,实测系统性低于非相干和(相干干涉损失),容差覆盖之 | ±12 dB |

**定标链注意**:dBm → W 为 $10^{(\mathrm{dBm}-30)/10}$(易错点:误写 +30 会
使信号弱 30 dB 全部没入噪声)。

## 6b. 下游处理链(radar_proc/)

从 echo cube 恢复每站视角的距离与角度,纯 FFT 路径。

### 6b.1 Range FFT(距离)

快时间维 FFT,距离轴(单站近似,双基差毫米级可忽略):

$$
r_k = \frac{k\,F_s}{N_{fast}}\cdot\frac{c}{2S}
$$

非相干积累所有 (slot, rx) 得距离谱 $P_{acc}(k)$,谱最大峰 + 抛物线内插给距离。

### 6b.2 TDM 重排与角度 FFT

slot 顺序 $s = p + n_{tx}(c_{ch}-1)$ 与 MATLAB 列主序 reshape 匹配:
$(k,\ s,\ q) \to (p,\ c_{ch},\ q)$。chirp 维平均(静止目标相位不变),tx/rx
两维加 hamming 后 2D FFT,fftshift 峰 + 抛物线内插得数字角频率。

**相位步进符号(实测验证)**:发射程与接收程的距离都随阵元坐标沿视线投影
缩短,两维步进**同为负号**(TX 竖臂沿 z 测俯仰、RX 横臂沿 y 测方位):

$$
\omega_{tx} = -k_0 d_{tx}\sin el,\qquad
\omega_{rx} = -k_0 d_{rx}\sin az \cos el,\qquad k_0 = \frac{2\pi}{\lambda} \tag{E17}
$$

峰位换算:

$$
\sin az \cos el = -\frac{\hat\omega_{tx}}{k_0 d_{tx}},\qquad
\sin el = -\frac{\hat\omega_{rx}}{k_0 d_{rx}},\qquad
el = \mathrm{asind}(\sin el),\ \ az = \mathrm{asind}\!\Big(\frac{\sin az\cos el}{\cos el}\Big)
$$

角度分辨率受虚拟孔径限制(TX 竖臂测俯仰 / RX 横臂测方位,波束宽几十度),补零
(`cfg.proc.nfft_az / nfft_el`,默认 256)细化谱线后配合抛物线内插,单目标
误差 < 0.2°。补零不改变物理分辨率(由孔径决定),只提升峰定位精度,并避免
3/4 点短序列直接内插的失真。**实现易错点**:chirp 平均后的数组须 squeeze
再 fft2,否则 RX 维落在第 3 维不被 fft2 变换,角度恒错。

### 6b.3 真值对照

每站对每目标中心投影(站系):$C^s = R^{\mathsf T}(C - p)^{\mathsf T}$,
$r = \|C^s\|$,$az = \mathrm{atan2}(y^s, x^s)$,$el = \mathrm{atan2}(z^s, \sqrt{x^s+y^s}^2)$。
main 打印估计 vs 真值表。

## 7. 文件结构

```
3DRadarCalibration/
├── config.m                    % 全部参数,唯一入口
├── main.m                      % 静态流水线 + 自检 + 存盘
├── main_motion.m               % 运动轨迹遍历入口(§9)
├── scene/
│   ├── Scene_radarPoses.m      % 6D 位姿 + L 型阵元坐标(E10/E16)
│   ├── Scene_targets.m         % 方体 6 面面元网格
│   ├── Scene_trajectory.m      % Z/C/M 轨迹弧长采样到帧网格(E18)
│   ├── Move_target.m           % 目标平移(§9)
│   └── Plot_scene.m            % 房间 + 方体面元 + 阵列放大插图
├── signal/
│   ├── Waveform_derive.m       % E1 E2
│   ├── Antenna_gain.m          % E12
│   ├── Echo_generate.m         % E3-E13 核心
│   ├── Check_echo.m            % 六项自检
│   └── Plot_antenna.m          % 单元方向增益诊断图
├── radar_proc/
│   ├── RangeFFT.m / AngleFFT.m / Process_cube.m   % 处理链(§6b)
│   ├── Truth_project.m         % 真值投影对照
│   ├── Det_to_points.m         % (r,az,el) -> 站系直角坐标(§10)
│   ├── Fuse_SVD.m              % 配对点云 SVD 刚体配准(§10)
│   ├── Plot_processing.m       % 每站距离谱 + 2D 角度谱
│   ├── Plot_tracking.m         % 轨迹遍历总图(§9)
│   └── Plot_fusion.m           % 点云配准叠加 + 残差图(§10)
└── results/                    % <date>_echo.mat / <date>_motion.mat / png
```

## 8. 不建模清单与风险

**不建模**:直耦泄漏、站间互扰、地面墙体多径、chirp 边缘效应(τ ~ ns ≪ Tc)、
帧内运动与多普勒相位斜率(§9 论证可忽略)、接收链幅相不一致。

| 风险 | 预案 |
|---|---|
| 镜面峰动态范围大(59 dBsm 量级) | 物理正确;非镜面面元由 cos²θ_b 压回 |
| 面元数放大(2 m 方体 → 9600) | 分块广播已控制临时量;必要时提高 facet_size |
| phase_mode='random' 下自检 5/6 失效 | 自检固定用 mirror 构造;random 仅目标行为差异 |

## 9. 运动目标(轨迹遍历,main_motion)

把轨迹看成**多个静态点的遍历**:外层循环逐帧把目标平移到轨迹点,每帧内部
与静态流程完全一致(Echo_generate → Process_cube → Truth_project)。循环
内部不画图;画图统一在循环后由 `cfg.sim.plot` 开关控制(静态 main 的出图
同受此开关守卫,默认 true 行为不变)。

### 9.1 帧周期与停-走冻结(E18)

一次 cube 采集的帧周期:

$$
T_{frame} = n_{tx}\,N_{chirps}\,T_c = 9.6\ \text{ms} \tag{E18}
$$

轨迹按帧网格弧长匀速采样(Scene_trajectory),第 $n$ 帧目标冻结在
$\boldsymbol{p}_n = \mathrm{traj}(t_n)$,回波按静态模型在该位置完整重建
(可见性、方向图、镜面闪等几何效应随位置自然演化)。近似量级:

| 位移 | $v = 1$ m/s | $v = 2$ m/s | 对比 |
|---|---|---|---|
| 帧间(逐帧重建保留) | 9.6 mm | 19.2 mm | ≪ 距离分辨率 10 cm |
| 帧内(冻结忽略) | 50 µm | 100 µm | ≪ 分辨率 |

### 9.2 多普勒可忽略论证(E19)

帧内运动的一阶效应是双基多普勒及 up-chirp 距离-多普勒耦合:

$$
f_D = \frac{(\boldsymbol{u}_{tx} + \boldsymbol{u}_{rx})\cdot \boldsymbol{v}}{\lambda},\qquad
\Delta R = \frac{c\, f_D}{2S} = \frac{v_r\, c}{\lambda\, S} \approx 2.6\ \text{mm} \times \frac{v_r}{\text{m/s}} \tag{E19}
$$

距离峰被推移 $\Delta R$:$v_r = 6$ m/s 时约 1.6 cm,仍远小于 10 cm 分辨率。
故不建模 chirp 维多普勒相位斜率,只测 (r, az, el),适用速度上限约
10 m/s($\Delta R$ 逼近分辨率量级)。

### 9.3 轨迹形状与实现

- 形状锚点(Scene_trajectory,房间系 xy,高度 = `cfg.target.targets(1).pos(3)` 联动):
  **Z** 三段折线 [2.5 1.5]→[5.5 1.5]→[2.5 4.5]→[5.5 4.5];**C** 圆弧
  (圆心 [4 3]、r = 1.5、弧 60°→300°,开口朝 +x);**M** 四段折线
  [1.8 1.8]→[2.9 4.2]→[4.0 1.8]→[5.1 4.2]→[6.2 1.8]
- 弧长匀速采样到帧网格,超出 `n_frames_cap` 截断;末端冻结在终点
- 平移用 Move_target(pos 与 facet.pos 同移,phase 不动——散射相位是
  目标属性;等效于改 cfg 里目标位置重跑 Scene_targets,但不重抽 random 相位)
- 检测判定:峰SNR($P_{acc}(k)/\mathrm{median}$)≥ `cfg.motion.snr_th`(12 dB)
- 自检不变:Check_echo 验证的静态模型在每帧内严格成立

## 10. 双站点云融合自校准(SVD,main_motion 内置)

### 10.1 站系点云累积(E20)

逐帧检测 $(r, az, el)$ 按站系球坐标转直角坐标并逐帧累积(Det_to_points,
与 Truth_project 投影互逆):

$$
\boldsymbol{P}_i(n) = r_n\begin{bmatrix}\cos el\cos az\\ \cos el\sin az\\ \sin el\end{bmatrix} \tag{E20}
$$

两站各得一条点云——同一物理轨迹在两站系的表达,帧序即配对关系(同一帧
两站观测同一目标位置),故为**已知对应**的配准,无需 ICP 最近邻迭代。
双站同时检出的帧(`both = all(detected,2)`)进入配准。

### 10.2 相对位姿真值(E21)

$C_{s1} = R_1^{\mathsf T}(C - p_1)$,$C_{s2} = R_2^{\mathsf T}(C - p_2)$,消去 $C$:

$$
C_{s1} = \underbrace{R_1^{\mathsf T} R_2}_{R_{21}}\, C_{s2} + \underbrace{R_1^{\mathsf T}(p_2 - p_1)}_{t_{21}} \tag{E21}
$$

融合目标即估计 $\hat R \approx R_1^{\mathsf T}R_2$、$\hat t \approx t_{21}$;
配合站 1 的先验位姿即得两站绝对外参(自校准)。

### 10.3 Kabsch/Umeyama SVD 解(E22)

配对点云去质心 $X = C_2 - \mu_2$,$Y = C_1 - \mu_1$,互协方差
$H = X^{\mathsf T}Y$,SVD $H = U\Sigma V^{\mathsf T}$:

$$
\hat R = V\,\mathrm{diag}\big(1,1,\mathrm{sign}\det(V U^{\mathsf T})\big)\,U^{\mathsf T},\qquad
\hat t = \mu_1 - \hat R\,\mu_2 \tag{E22}
$$

轨迹点云位于 $z$ 恒定平面(Z/C/M 形状非旋转对称),平面法向约束倾角、
面内形状约束转角,6 自由度均可观。配准用 3σ 两遍法:首遍全配对解,
残差 $> 3\times\mathrm{RMSE}$ 的配对剔除后重解(实测剔除量少,精度由逐点
量测噪声主导)。实测精度(534–551 配对帧,单点横向噪声 ~2.6 cm @
0.3°/5 m):旋转误差 0.28–0.45°,平移误差 68–84 mm。

### 10.4 构造房间系显示(Roomframe_plane,Plot_fusion)

站 1 绝对位姿不可观(全局 6 自由度规范自由)。目标同水平面运动先验下,对站 1
内点云做平面拟合,数据可观 3 个自由度:法向定房间竖直方向(2),偏距定雷达离
目标平面的垂直距离 d(1);其余 3 个为纯约定——偏航取站 1 波束的水平投影方向
(房间 x 轴),原点取站 1 正下方,目标平面取 z = 0。组装 $R_1$ = $[\hat x\ \hat y\ \hat z]^{\mathsf T}$
(行 = 房间轴在站 1 系表达),$p_1$ = $[0,0,d]$,$C_w = R_1 C_{s1} + p_1$。
该系仅用于显示(Plot_fusion),不是墙对齐的物理房间系;配准与绝对外参推导
(式 E21)不受影响。当前默认配置下 d = 站 1 高度 − 目标平面高度 = 2.7 m。
