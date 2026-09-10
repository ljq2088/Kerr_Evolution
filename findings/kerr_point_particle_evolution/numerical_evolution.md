# Kerr 点粒子时域演化

- 来源任务：T007 / W02.3
- 状态：owner accepted
- 连续理论：[主推导](theory/DERIVATION.md)
- 正式 solver：[kerr_point_particle_solver](../../code/kerr_point_particle_solver/README.md)
- 公共工具：[kerr_waveform_tools](../../code/kerr_waveform_tools/README.md)
- 正式数据：[SXS:BBH:0305 medium bundle](../../data/kerr_point_particle_evolution/sxs0305_medium/manifest.json)
- 独立证据：[最低分辨率跨代码对比](evidence/cross_code_low_resolution.md)

## 数值问题

对固定 $m$，代码在

$$
(T,R,y),\qquad R=\frac{M^2}{r},\qquad y=-\cos\theta
$$

上演化 peeling variable $\psi_{4,m}$ 的 two-field method-of-lines 系统

$$
P=A\psi_{,T}+B\psi_{,R}+D\psi,
\qquad
\dot U=\mathcal F(T,U),\qquad U=(P,\psi).
$$

空间使用七点差分及 horizon、$\mathscr I^+$、axis 的实际闭合行；时间使用 classical RK4。每个 RK stage 在对应时间重新采样唯一 worldline 并评价四块 source。点粒子穿过 horizon 后 source 关闭，trajectory 本身继续保留到 terminal event。

默认 production 设置为

$$
(N_R,N_y)=(512,129),\qquad
\sigma_R=4\Delta R,\qquad
\sigma_y=4\Delta y,\qquad
\tau_{\rm on}=20M,
$$

并在 $\mathscr I^+$ 投影得到 complex $\psi_{4,\ell m}(T)$。内部 RK 步长由 wave/source 条件保守选择，和 $\Delta T_{\rm out}$ 明确分离。

## 数值接受依据

在 $\chi=0.8$、$(\ell,m)=(2,2)$、$(N_R,N_y)=(128,49)$ 上，当前实现与已经完成 convergence/QNM 检查的独立目标代码直接比较，没有时间、相位或振幅对齐。在 $107M\le T\le259.6M$：

$$
\epsilon_{\rm cross}=0.454\%,
\qquad
\|\Delta\psi_{4,22}\|_\infty=1.49\times10^{-4},
\qquad
|\mathcal O|=0.9999907.
$$

目标代码自身最低档到下一档的相对 RMS 差异为 $5.76\%$，约为跨代码差异的 $12.7$ 倍。Owner 据此接受当前 solver，并取消原计划中其余数值 validation gates；被取消的检查保留在历史 plan，不改写为“已通过”。

## SXS:BBH:0305 输出

代表性背景为

$$
\frac{M_f}{M}=0.952032939704,
\qquad
\chi_f=0.6920851868180025.
$$

$(2,2)$ 与 $(4,4)$ 分别完成 medium 演化，均含 3587 个 $\Delta T=0.1M$ 样本，终点 $T_1=358.7M$。正式 bundle 同时保留：

- 原始 complex $\psi_{4,22}$ 与 $\psi_{4,44}$；
- 完整 FFI $H_{22}=rh_{22}$ 与 $H_{44}=rh_{44}$；
- 每个 mode 的实部/绝对值基础图；
- SXS h22/h44 输入、精确自旋 ISCO flux、plunge 初值、run metadata 和全部 hashes。

原始 $\psi_4$ 是权威 solver 输出；strain 是受控派生产品，不能替代原始 mode。

## Strain 与边界

C2 使用完整实际时间序列和 P19 fixed-frequency integration：

$$
\widetilde H_{\ell m}
=-L^2\frac{\widetilde\psi_{4,\ell m}}
{\max(|\omega|,\omega_0)^2},
\qquad
\omega_0=0.75|m\dot\Phi_p(T_{\rm stable})|.
$$

输出保留 transform、analysis、recommended intervals 和实际 taper。现有 A1 strain sign/normalization bridge 仍按 artifact metadata 标记为 `conditional-open`；owner 接受 W02 数值结果不等于一次未执行的独立 convention proof。Sampling 状态记录为 `owner-waived`。

## 下游接口

W03/W04 直接读取正式 complex $\psi_4/H$ 和事件、质量、自旋、convention/provenance；W05/W06 复用公共 Kerr background、ISCO/transition、worldline、time jets、harmonic 与 artifact 接口。QNM filter、近视界 radial cache 和 SPA 仍由各自任务建立，不由 W02 迁移预先定义。
