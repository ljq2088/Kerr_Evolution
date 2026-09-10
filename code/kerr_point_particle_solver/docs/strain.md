# C2 FFI strain 实现

- 状态：owner accepted
- 输入：完整、均匀采样且 `status=complete` 的 complex $\psi_{4,\ell m}(T)$
- 输出：$H_{\ell m}=r h_{\ell m}$；不除以 luminosity distance，也不恢复物理尺度
- Sign/normalization：`conditional-open`；C2 不独立关闭 A1 strain bridge
- 正式实例：[SXS:BBH:0305 medium bundle](../../../data/kerr_point_particle_evolution/sxs0305_medium/README.md)

## 变换

实现采用 P19 fixed-frequency integration：

$$
\widetilde H_{\ell m}^{\rm FFI}(\omega)=-L^2\frac{\widetilde\psi_{4,\ell m}(\omega)}{\max(|\omega|,\omega_0)^2},
\qquad \omega_0=0.75|m\dot\Phi_p(T_a)|.
$$

FFT 使用 mode 文件的完整实际时间序列。$T_a=T_{\rm stable}$ 是分析起点而不是
截断点，$T_b=T_1$；两端采用 $10M$ quintic taper。输出保存完整 $T$、complex
$H_{\ell m}$、实际 taper、transform/analysis/recommended intervals、cutoff、
输入 hash、source events、背景、环境与原始 mode metadata。

## Workflow 合同

`save_strain_lm=True` 会在完整 PDE run flush 后自动执行 C2。四种输出请求均有
测试：只保存 $\psi_4$、只保存 $H$、两者都保存或两者都不保存。只请求 strain 时，
中间 $\psi_4$ 在 staging 中支持恢复；只有 $H$ 写入并经 loader 复核后才删除。
短 test run 不允许生成标为 `complete` 的 strain。

新输出采用 `c2-H-v2`，reader 兼容历史 `c2-H-v1`。Sampling 状态按 owner 决定
记录为 `owner-waived`；cutoff variation 是 systematic diagnostic，不反向改变
baseline。解析单频、低频污染、完整早期信号、保留频带 forward reconstruction、
非法输入、原子写入、loader 和 provenance 均由自动测试覆盖。
