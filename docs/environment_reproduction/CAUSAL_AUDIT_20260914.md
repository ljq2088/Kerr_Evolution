# Kerr 偏差的因果核查：本轮结果

2026-09-14。WSL 项目 `/home/ljq/code/kerr-hyperboloidal`，分支 `codex/reproduce-kerr-environment`。

**结论：尚未定位到能解释并消除 32.6552% 总视界通量差异的根因。** 本轮找到了并修复了一个诊断原型错误；同时完成了有源匹配与独立径向边界求解，量化否定了此前最优先的低阶规范幅度假设。不能把这两件事混称为“修复 Kerr 通量”。

## 1. 确实找到的代码错误及其范围

`paper_jump_basis.angular_jet` 原来从 `separated_jet` 提取所谓角函数。后者返回的是带 Kinnersley 因子的分离场：

\[
F_s=R_s S_s/\zeta^{|s|-s}.
\]

对于负自旋，取其角向 Taylor 系数并不等于纯角函数 \(S_s\)。这使上轮新增的局部匹配矩阵带入了错误的径向因子。现在直接用角函数值、一阶导数和角向 ODE 构造纯角 Jet，径向导数严格为零；还增加了角向 Teukolsky–Starobinskii 恒等式检验。

旧测试只比较两个共享错误角函数的计算路径，因此没有发现这个错误。新增测试直接对照角本征函数，避免同源错误相互验证。

**这个函数只属于新建的匹配诊断原型，未被历史环境源和通量调用。** 修复它不会改变此前 Kerr 生产结果，也不能解释 Schwarzschild 与 Kerr 的区别。底层 SWSH 本征函数本身不是这次查到的错误。

上轮条件数结论更正如下（问题参数 a=0.8771530275949366、r0=20M）：

| 局部系统 | 原错误值 | 更正值 |
|---|---:|---:|
| m=1、L=1 | 52.9834 | 5040.4188 |
| m=1、L=4 | 53.0260 | 5041.2500 |
| m=2、L=2 | 93.7281 | 374.3343 |

偶极缺两条常规约束、补入 rho h_l+m+ 后恢复满秩的结构结论仍成立；这不是旧生产算法漏条件的证据。

## 2. 独立有源匹配：从物理粒子源求规范幅度

新增 `paper_sourced_matching.py`。右端包含已知 Weyl 跳跃、trace 跳跃、对角 kappa 的有源高阶导数及非对角 kappa；未知量为每阶的

\[
x_\ell=([P_{-1}],[P_{-1}'],[\kappa_{\ell\ell}],[\kappa_{\ell\ell}']) .
\]

把三个显式标架分量投影到球谐上，施加场连续与粒子应力张量给定的导数跳跃。常规选 h_l+l+、h_m+m+；m=j=1 用 rho h_l+m+ 替代不存在的 spin-2 球谐条件。行列缩放后用带列主元的复 QR 求解，不从环境通量反推系数。

在公开表格参数 r0=6、a=0.6、m=2、L=10 下得到：

| ell=2 的量 | 独立源匹配结果 |
|---|---:|
| [P_-1] | -312.0748083698 + 93.5611288136 i |
| [P_-1'] | 66.5946478109 - 27.2132034014 i |
| [kappa] | -27.6725775375（虚部约 9e-11） |
| [kappa'] | 35.0620909506（虚部约 -3e-10） |

低阶结果与公开表格打印精度相容；Schwarzschild 表格另外纳入自动回归。这里没有宣称整个 L=10 截断都已验收：最高保留阶的未使用方程残差明显，不能由被求解方程的微小残差推断全系统精度。

在实际问题 a=0.8771530275949366、r0=20、m=1 下，以新匹配结果替换旧偶极跳跃，并将差值恢复为满足 IN/UP 边界的齐次径向修正。旧 spin-1 向量转换到论文 P_-1 基底后的相对拟合残差为 7.94e-15。该拟合仅作基底换算，没有拟合物理通量。

将规范向量修正变成标量特解

\[
\delta f=-\delta\xi^{a*}\partial_a\Phi_c,
\]

其中只共轭度规模态，云不共轭。再用分段 Green 恒等式计算视界复振幅修正，保留轨道两侧项和通量干涉：

\[
\delta Z_H=\left[\frac{\Delta(U\delta f'-U'\delta f)}{W}\right]_{r_{in}}^{r_{out}}
+B(r_0^-)-B(r_0^+),
\qquad
\frac{\delta F_H}{F_H}=\left|1+\frac{\delta Z_H}{Z_H}\right|^2-1.
\]

| 新匹配内部阶数 | delta Z_H | 相对通量变化 |
|---|---|---:|
| L=4 | 2.33371195e-8 + 1.00368336e-8 i | -4.55782662e-7 |
| L=6 | -3.54877153e-8 - 3.11363723e-8 i | +3.99119352e-7 |

原偶极 Z_H=-0.0500392802758+0.0337036736095i。上述变化换成百分比约为 -0.0000456% 和 +0.0000399%，与总视界偏差 32.6552% 相差约六个数量级。这里传播的是 ell=1 齐次规范修正，不是新全模态总通量。

低频高阶导数存在精度损失；扩展 Jet 系数并没有使双精度角函数、径向输入和 QR 变成任意精度计算。L4/L6 小修正符号不同，不作为严格收敛误差界。尽管如此，现有直接计算不支持“旧偶极跳跃错了约百分之十，从而造成主通量差异”的解释。

数据：`paper_sourced_jumps_*.json`、`paper_dipole_closure.json`。

## 3. 独立检查 kappa 的全域特解与两端边界

仅比较跳跃还可能漏掉特解/齐次解分拆。本轮进一步新增 `paper_kappa_radial.py`，不调用质量导数法，直接解 2023 方程：

\[
(\Delta h')'+V_\ell h=0,\qquad
(\Delta\kappa_{\ell\ell}')'+V_\ell\kappa_{\ell\ell}
=\tfrac12(r^2+a^2\gamma_{\ell\ell})h.
\]

无穷远取出射特解渐近展开，视界取入射 Frobenius 展开，分别积分耦合的 h、kappa ODE。端点特解的自由齐次系数先置零，再由给定的粒子处跳跃求齐次补偿。该任意中间选择会被补偿吸收，不根据通量决定。

与既有

\[
\kappa_{23}=\tfrac12\partial_{\mu^2}h-\frac{i}{2\omega}\chi_c
\]

的对角径向分量比较。使用同一物理 trace 源和规定跳跃，比较场值及一阶径向导数。

| 背景/模态 | 新方法远端 | 检查点上的最大相对差 |
|---|---:|---:|
| a=0、r0=6、ell=m=1 | 4000M | 1.85e-10 |
| a=0.6、r0=6、ell=m=2 | 4000M | 4.55e-12 |
| a=0.877153、r0=20、ell=m=1 | 2000M | 1.673e-11 |
| 同上 | 4000M | 1.673e-11 |
| 同上 | 8000M | 1.673e-11 |

Kerr 检查覆盖 r/M=1.481、3、10、19.9、20.1、40、200。结果支持两种方法给出同一个对角特解加齐次补偿，而非仅仅满足同一个真空方程。它不是作者原始径向数组的比较，也不是全域严格误差界。

数据：`paper_kappa_radial_audit.json`。

## 4. spin-2 与参考图的补充核验

在 Schwarzschild、a=0.6 及问题 Kerr 参数下，将 2023 显式 spin-2 公式与现有协变重构的三个分量直接比较。共享 Weyl 径向输入，重构公式独立。五个点的最大相对差为 1.51e-9，其余为 2e-13 至 2.34e-11。该检查不覆盖所有十个分量及全域高阶精度。

数据：`paper_tensor_comparison.json`。

再次查看 arXiv v1 图2与期刊作者接受稿的同图，视觉上未发现曲线替换。没有据此声称两个 PDF 向量路径严格一致，也没有发现可直接解释偏差的勘误。视界通量的 4Mr+ 因子与论文的模态公式相符；不能修改面积因子来追随曲线。

## 5. 现在能够与不能够说什么

- 可以说：新局部匹配原型有明确的角函数包装错误，已修复并补上独立回归。
- 可以说：此前最优先怀疑的低阶规范跳跃差异，传播到通量后量级不足；kappa 特解边界的独立求解也与现有结果吻合。
- 不能说：该原型错误造成了历史 Kerr 偏差，或者本轮已修复 32.6552%。生产通量仍为原值。
- 不能说：Schwarzschild 某个点接近就证明整个源/重构模块对所有 Kerr 参数正确；也不能因剩余偏差就断言作者代码有错。

下一次能真正区分原因的比较对象应当是作者在同一 Kerr 参数下的**已归一化云径向函数、m_g=1 度规分量、scalar00 投影源和 Z_H 复振幅**，或一条不共享当前源构造的完整实现。仅有最终曲线，无法直接区分输入场归一化、源组装、输出约定或作者数值处理。当前尚无这些作者中间量。高阶重构污染也仍待独立修复，但已有偶极归因显示它不能被直接认定为主偏差。

## 6. 复现与验证

在项目根目录设 `PYTHONPATH=src OPENBLAS_NUM_THREADS=1`，使用 `.venv/bin/python`：

```bash
.venv/bin/python src/report_paper_sourced_matching.py --a 0.8771530275949366 --r0 20 --m 1 --ellmax 4 --quadrature 20 --extended --output outputs/paper_matching_check.json
.venv/bin/python src/report_paper_dipole_closure.py
.venv/bin/python src/report_paper_kappa_radial.py
.venv/bin/python src/report_paper_tensor_comparison.py
.venv/bin/python -m pytest -q tests/test_paper_jump_basis.py tests/test_nonstatic_paper_projection.py tests/test_static_tetrad_projection.py
```

相关回归 8 项通过。匹配诊断的高阶未使用残差与精度限制保留在结果中，没有以测试通过替代论文一致性验收。

## 原始资料

- [环境论文 arXiv:2501.09806v1](https://arxiv.org/abs/2501.09806v1)
- [作者接受稿，期刊论文](https://eprints.whiterose.ac.uk/id/eprint/227286/8/Self_force_and_environments%20(2).pdf)
- [2023 Kerr Lorenz 重构](https://arxiv.org/html/2306.16459v3)
- [2024 直接 Lorenz 重构](https://arxiv.org/html/2406.12510v3)
