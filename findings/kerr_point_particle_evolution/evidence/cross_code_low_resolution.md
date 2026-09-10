# 最低分辨率跨代码波形对比

- 状态：`owner_accepted`
- 物理参数：$M=1$、$\chi=0.8$、$(\ell,m)=(2,2)$
- 空间分辨率：$(N_R,N_y)=(128,49)$
- 对齐：共同绝对 $T$；没有时间、相位或振幅对齐
- 目标结果：[独立代码最低分辨率波形](/Users/indigo/codes/Direct_Wave_in_Higher_order/Codes/NewPointParticlesm2Evolution/runs/sminus2_low_chi0.8_m2_20260810_183732/psi4_l2_m2_R0.npz)

## 结果

目标仓库的收敛分析使用 $107M\le T\le259.6M$。在同一窗口内，

$$
\frac{\operatorname{RMS}
\left(\psi_{4,22}^{\rm current}-\psi_{4,22}^{\rm target}\right)}
{\operatorname{RMS}\left(\psi_{4,22}^{\rm target}\right)}
=4.54\times10^{-3},
$$

$$
\left\|\psi_{4,22}^{\rm current}
-\psi_{4,22}^{\rm target}\right\|_\infty
=1.49\times10^{-4},
\qquad
\left|\mathcal O_{\rm complex}\right|=0.9999907.
$$

两套代码都在 $T=207M$ 达到峰值。峰值分别为

$$
\left|\psi_{4,22}^{\rm current}\right|_{\rm peak}
=0.0225306931,
\qquad
\left|\psi_{4,22}^{\rm target}\right|_{\rm peak}
=0.0225305564.
$$

目标代码自身从最低分辨率 $(128,49)$ 到下一档 $(256,65)$ 的同窗差异为

$$
\epsilon_{\rm res}^{\rm target}=5.76\%,
\qquad
\|\Delta\psi_{4,22}\|_\infty=1.49\times10^{-3}.
$$

当前跨代码差异分别小约 $12.7$ 倍和 $10.0$ 倍，明显低于目标实现自身的该档离散差异。按“两个独立实现应在既有最低分辨率误差内重合”的判据，本次对比通过。

## 约定差异

目标旧 run 没有 startup ramp；当前代码用 $\tau_{\rm on}=20M$ 的 quintic ramp。因此早期 transient 明显不同，不能把 $T<107M$ 的全时段 RMS 用作代码主体是否一致的判据。目标旧 run 还使用 $\Delta T_{\rm RK}=0.01M$，并让 Gaussian 尾在粒子过视界后保留到 $T=202.969M$；当前代码使用 $\Delta T_{\rm RK}=0.007142857M$，在 $T=198.196M$ 关闭源。这些差异保留在比较中，没有进行补偿。

从 $107M$ 到 source crossing 前的相对复数 RMS 为 $0.190\%$；ringdown 段为 $0.712\%$。因此上述约定差异没有破坏主波形与晚时波形的一致性。

## 结论

目标实现已经具有既有 convergence 与 QNM 检查；当前实现由另一套代码路径独立构造。二者在相同物理参数和最低空间分辨率下的主波形远比目标实现自身的分辨率误差更接近，峰值时间与振幅也一致。这是当前新代码数值实现成功的强证据。

Owner 依据这一结果接受当前 solver，并决定不重复执行原计划中的其余数值 validation gates。按 T007 迁移决定，过程性 comparison run、checkpoint、JSON 和图片在关键数值写入本记录后删除；目标仓库原始参考路径保留在上方。
