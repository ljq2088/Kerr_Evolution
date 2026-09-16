# scalar00 轨道拼接与 Lorenz 接触源审计（2026-09-16）

**实际测量不支持轨道漏掉的 Lorenz 接触项解释现有视界偏差。** 在 \(a/M=0.8771530275949366\)、\(r_p/M=20\)、\(\alpha=0.3\)、\(m_g=-1\) 下，直接计算生产完整度规的轨道两侧极限、做 \((\ell,m)=(0,0)\) 角向投影，并将可能的接触源传播到视界。五组控制给出 \(|\Delta Z_H/Z_H|<1.8\times10^{-6}\)；其中直接计算负 \(m_g\) 的生产路线为 \(5.20\times10^{-7}\)。这比需要解释的约 12% 幅度差小近五个数量级。

这不是全度规匹配的无条件通过证书。有限 \(L=18\) 的部分逐点 BL 度规跳跃仍明显，且高阶余项没有单调收敛到零；被排除的是它们经实际目标收缩、投影后在 scalar00 中形成足够大轨道接触源的解释。本次未修改生产模块、未把诊断接触项加入物理结果，也没有重跑作者 Mathematica 程序。

## 计算对象和符号

在 \((\Box-\mu^2)\Phi_0=0\) 时，完整一阶标量方程是

\[
 (\Box-\mu^2)\delta\Phi=h^{ab}\nabla_a\nabla_b\Phi_0
 +C^b\nabla_b\Phi_0,\qquad C^b=\nabla_a\bar h^{ab}.
\]

对于直接按内、外径向解拼成的度规，

\[
 C^b=C_-^b\Theta(r_p-r)+C_+^b\Theta(r-r_p)
 +[\bar h^{rb}]\delta(r-r_p),\qquad [h]=h_+-h_-.
\]

因此原有 \(h:\nabla\nabla\Phi_0\) 管线可能遗漏的轨道项，只需计算

\[
 J_\delta^{00}=2\pi\int_{-1}^1 S_{00}(x)\,
 \Sigma(r_p,x)[\bar h^{rb}](r_p,x)\nabla_b\Phi_0(r_p,x)\,dx,
\]
\[
 \Delta Z_H={R_{\mathrm{Up}}(r_p)\over W}J_\delta^{00},\qquad
 W=\Delta(R_{\mathrm{In}}R'_{\mathrm{Up}}-R_{\mathrm{Up}}R'_{\mathrm{In}}).
\]

此处已经使用了径向方程 \(\partial_r(\Delta\partial_rR)+VR=J\) 的源约定，不能再额外乘除 \(\Delta\)。度规在轨道处先取协变分量极限，再用同一个背景逆度规升指标并做 trace reversal。所有 \(rt,r\theta,r\varphi,rr\) 收缩均保留；没有只检查 \(h_{rr}\)。

生产中所有非静态扇区都包括在内：spin 2、完整双手性 spin 1、trace、compact chi 和 kappa，逐个 \(\ell_g=1,\ldots,18\) 累加。没有删 vector。第一路线由正 \(m_g\) 全度规共轭得到负 \(m_g\)，云保持原来 \(m_b=1\) 的复场；独立控制则直接调用 `nonstatic_metric(...,m=-1)`，不共轭输出，符合基础强迫源生成器的直接路线。

## 外推、投影和实际结果

在每个 Gauss–Legendre 角点，分别取 \(r=r_p\pm\epsilon\)，用 order-8 Taylor jets 提供度规及其前两阶径向导数，外推

\[
 h(r_p^\pm)=h(r_p\pm\epsilon)\mp\epsilon h'(r_p\pm\epsilon)
             +\tfrac12\epsilon^2h''(r_p\pm\epsilon).
\]

分别改变角向点数、将 \(\epsilon\) 减半、将 Taylor 系数改用 `clongdouble`，并改用直接负 \(m_g\) 路线。extended 仅提高本地 Taylor 运算精度；输入的径向/角向解仍是原后端精度，不能称作全流程任意精度计算。各组均用相同云归一化和基准 Green 定义。

| 度规计算路线 | 角点数 | \(\epsilon/M\) | \(|\Delta Z_H/Z_H|\) | 假设加入接触项的相对通量变化 |
|---|---:|---:|---:|---:|
| 共轭正 m，double | 12 | 1.0e-03 | 1.7671e-06 | -3.9491e-07 |
| 共轭正 m，double | 24 | 1.0e-03 | 9.3230e-07 | 1.5558e-06 |
| 共轭正 m，double | 24 | 5.0e-04 | 3.0497e-07 | 5.3156e-07 |
| 共轭正 m，extended Taylor | 24 | 1.0e-03 | 2.8520e-07 | 4.6537e-07 |
| 直接负 m，double | 24 | 1.0e-03 | 5.2020e-07 | 9.6619e-07 |

最后一行的直接负 \(m_g\) 结果具体为

\[
 J_\delta=(-1.94085+1.16600i)\times10^{-7},\qquad
 \Delta Z_H=(-1.75017+2.57585i)\times10^{-8},
\]
\[
 Z_H^{\mathrm{base}}=-0.0496101415924+0.0335053762021i.
\]

通量变化按 \(|1+\Delta Z_H/Z_H|^2-1\) 计算，仅用于量化假设遗漏项的影响，没有作为物理修正发布。

目前完整 scalar00 通量为 \(-6.70075034\times10^{-5}\)，论文图读**总**视界通量为 \(-5.17637534\times10^{-5}\)。即使暂不计其余同为负的 Kerr 模态，仅把 scalar00 幅度缩至这个总量也至少需要

\[
 {|\Delta Z_H|\over|Z_H|}\ge
 1-\sqrt{|F^H_{\mathrm{paper,total}}|/|F^H_{00}|}
 =0.12107633.
\]

故这里不存在“微小接触项因为相消而产生 12% 幅度变化”的数值迹象：已经直接比较最终复振幅响应。

## 收敛限度：不能把小投影说成所有分量连续

以 q24 的 extended 控制为例，逐项累加到 \(L=1,2,3,4,8\) 的 \(|\Delta Z_H/Z_H|\) 分别约为 \(5.40\times10^{-5}\)、\(8.47\times10^{-6}\)、\(7.62\times10^{-8}\)、\(4.26\times10^{-9}\)、\(7.47\times10^{-9}\)。低阶残差在补足耦合后迅速变小；更高阶又出现数值/截断余项，不能把最终 \(10^{-7}\) 的具体末位当作收敛的非零物理值。

L18 点值度规跳跃的最大 BL 分量 Frobenius 范数约为 0.10–0.23（各组角点不同），不是零，且不随这些控制单调下降。该量依赖坐标分量及角点，不能当规范不变量；它也说明不能用本实验宣称整张轨道球面上的所有度规分量已经验证连续。scalar00 所需的是上面特定的 trace-reversed 法向跳跃与云梯度收缩后的角向投影；这两类判据不能混用。

轨道两侧的光滑 Lorenz 残差也保存了。q24 double 的 \(\Sigma C^b\nabla_b\Phi_0\) scalar00 投影约 \(10^{-8}\)；extended Taylor 后两侧分别约 \(5.14\times10^{-11}\)、\(3.28\times10^{-11}\)，支持浮点消减误差而非大的光滑漏源。它们只是轨道附近的点值，**不是**全径向区间的积分误差界；全域光滑算子的独立控制见 `direct_kg_variation_audit_20260916.json`。

## 已有 full_gauge_closure 能证明的范围

`full_gauge_closure_L8_precise.json` 将重新求得的有限 L gauge jumps 与旧值之差传播到 scalar00，给出相对通量变化 \(-5.05\times10^{-12}\)；带六阶截断保护的 `trusted_gauge_closure_L10.json` 为 \(3.18\times10^{-11}\)。这些结果约束的是 **两套 gauge 匹配系数之间的差**，共享曲率与 trace 输入，不能单独证明原完整度规的 \([\bar h^{rb}]\) 为零。本次直接单侧度规检验补充了这一缺口。

在当前目标模态和现有分辨率下，轨道局部接触项无法解释原 32.7% 总通量偏差。后续不应靠加这个微小项或删去 vector 去拟合论文；剩余追查仍需回到尚未独立闭合的输入、背景/单位约定或完整源与原始论文数据之间的差异。

## 复现文件

- `src/report_orbit_source_matching_contact.py`：逐侧全度规、Taylor 外推、trace reversal、scalar00 投影和 Green 响应；包含输入与依赖哈希。
- `src/report_orbit_source_matching_direct_negative.py`：独立直接负频率/负 m 控制，复用保存和投影部分；仅替换诊断 worker，不修改任何物理函数。
- `src/report_orbit_source_matching_summary.py`：只读汇总五组原始记录，包含局部 Lorenz 收缩。
- `docs/environment_reproduction/orbit_source_matching_summary_20260916.json`：五组结果、输入 SHA-256、逐 L 收敛信息；原始文件名以 `orbit_source_matching_contact_` 开头。

基准 `fresh_20260915_L18_mg-1_sl0.json` 在本次运行时通过当前生产源的 provenance 校验。上述测量针对该基准的 unit cloud mass、scalar00、精确阈值 Kerr；不推广到其他 m、Schwarzschild 或全空间规范证明。
