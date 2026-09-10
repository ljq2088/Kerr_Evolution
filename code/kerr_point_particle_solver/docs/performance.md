# RK4 演化性能优化记录

- 状态：owner accepted
- 范围：fixed-$m$ source/RHS 的等价代码优化
- 基准：$\chi=0.8$、$m=2$，相同 Python/NumPy/SciPy 环境；计时不包含 system assembly
- 计划：[RK4 演化性能优化计划](../../../workspace/plan/W02_point_particle_evolution/numerical/appendix/performance_optimization.md)

## 1. Profiling 结论

原始 smoke 实现中，单次 source-on RHS 约为 $3.87\,\mathrm{ms}$，其中 source 为 $3.70\,\mathrm{ms}$，占 $95.8\%$；纯场 RHS 约为 $0.25\,\mathrm{ms}$。完整 run 在 source active 时约需 $2.2\,\mathrm{s}/M$，而 source off 后达到约 $7.2M/\mathrm{s}$。因此本轮首先处理 source，而没有改变 $\Delta T_{\rm RK}$、网格、Gaussian、RK4 或输出采样。

## 2. 已采用的实现

P1 缓存了只依赖 $(R,y,M,a,L,m)$ 的 source 几何因子；同一 stage 的三个 tetrad amplitudes 共用一次 Gaussian 与其 time-jet 因子；source off 后复用只读零源。视界正则总源仍按

$$
S_m^{\rm num}\big|_{R_H}
=-16\pi R\widehat\Sigma
\left(H_1+H_2+H_3+H_4\right)\big|_{R_H}
$$

组合，但正式路径复用当次已经得到的 $H_A$，不再为一条视界行重算整张网格。独立的视界 evaluator 保留在测试中。

P3 根据

$$
\frac{(R-R_p)^2}{\sigma_R^2}
+\frac{y^2}{\sigma_y^2}
\le 2\log\epsilon_g^{-1}
$$

确定截断 Gaussian 的非零区域，并让每次嵌套差分依据实际七点矩阵逐层扩展输出支撑。选择条件使用包含两层 stencil halo 后的局部面积与全网格面积之比：smoke 的局部 bookkeeping 得不偿失，因此使用 P1 全网格路径；low/medium/high 使用局部路径。

## 3. 计时结果

下表为五组重复计时的中位数；不同命令之间存在少量系统波动，数量级和相对结论稳定。

| $(N_R,N_y)$ | P1 full-grid RHS | 最终 hybrid RHS | 最终 source | 路径 |
|---|---:|---:|---:|---|
| $(128,33)$ | $1.00\,\mathrm{ms}$ | $0.98\,\mathrm{ms}$ | $0.81\,\mathrm{ms}$ | full-grid |
| $(256,65)$ | $3.15\,\mathrm{ms}$ | $2.14\,\mathrm{ms}$ | $1.62\,\mathrm{ms}$ | local |
| $(512,129)$ | $13.83\,\mathrm{ms}$ | $4.64\,\mathrm{ms}$ | $2.58\,\mathrm{ms}$ | local |
| $(1024,257)$ | 未单独测 P1 | $18.19\,\mathrm{ms}$ | $5.75\,\mathrm{ms}$ | local |

最终 code-hash 下的完整 smoke run 用时约 $2\,\mathrm{min}\,01\,\mathrm{s}$，原始 run 为约 $7\,\mathrm{min}\,32\,\mathrm{s}$，端到端墙钟提升约 $3.74$ 倍。最终 run 的 code hash 为 `4644823852dca147ce09e1e757dd911d73411b0a689db6788ffb39656fd7dbd5`，与当前源码一致。

## 4. 数值等价性

- P1 的多个代表时刻四块 source、总源、RHS 与 20 个 RK4 steps 均与优化前快照逐 bit 相同。
- 最终 code-hash 下的 1,061 点完整复数 $\psi_{4,22}(T)$ 与原始 smoke 输出逐 bit 相同，最大绝对 residual 为零。
- Local path 与 full-grid path 已在 $\chi=0.2,m=2$、$\chi=0.95,m=4$、多个分辨率和视界前时刻比较四个 blocks 与总源。
- Local sparse restriction 会改变浮点累加分组；已测整体相对差异约为 $10^{-13}$，测试门槛为 $2\times10^{-12}$。Medium 网格从同一非零状态执行 20 个 RK4 steps 后仍通过该门槛。
- 当前测试结果为 27 passed。

这些比较说明候选优化保持当前离散实现；它们不证明 continuum convergence、source normalization 或物理 waveform 已通过。

## 5. 未采用与后续边界

直接 NumPy stencil 比现有 CSR 慢约 $4$--$6$ 倍。Numba 单算符在 smoke/low/medium 没有稳定优势；融合四个场导数在 high 网格可把导数部分缩短约 $21\%$，但对完整 source-on RHS 的预计收益仅约 $7\%$，本轮没有为此新增 JIT 依赖。

P4 的完整 RHS 编译与内存峰值测量未执行。原计划中的额外 convergence 检查随后由
owner 取消；本节只支持已执行的等价优化与计时结论。
