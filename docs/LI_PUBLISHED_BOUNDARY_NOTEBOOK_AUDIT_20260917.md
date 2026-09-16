# Li 公开边界 notebook 审计

日期：2026-09-17。只读原始资料，新增隔离参考及诊断产物；未修改生产径向算法。本次没有重跑源积分。

## 确定结果与范围

作者仓库已取得并固定到 **`da78d2c96af5854e0feea0421af1ebb6c45c6bc0`**，提交时间 `2025-06-24T23:53:05-05:00`。其唯一受版本控制文件是 `Boundary Coefficients.nb`，165979 字节，SHA256 `706c5aa8a01e6e1bec123d3101819211cdf46c1e8855e3bfb06aaee9509ef24e`。此次 clone 的 `main` 及全部远端历史只有一个 commit。

该 notebook 只公开 **A1–A4 与 B1–B4** 的符号表达式。它没有原图外边界半径、近视界起点、积分精度、实际轨道参数或作图驱动程序。因此，不能由此认定 Li 的原图用了 `Rout=4000`，即使某个有限外边界控制实验碰巧更接近原图。

主来源：[作者仓库固定版本](https://github.com/dongjun826/EMRI-in-scalar-clouds/tree/da78d2c96af5854e0feea0421af1ebb6c45c6bc0)，[固定 notebook](https://github.com/dongjun826/EMRI-in-scalar-clouds/blob/da78d2c96af5854e0feea0421af1ebb6c45c6bc0/Boundary%20Coefficients.nb)，[Li v2 正文](https://arxiv.org/html/2507.02045v2)。本机论文 `reference.bib:6005–6007` 的 `boundary_coeffs` 直接指向这个仓库，故是论文明确引用的补充产物。

## 精确前因子和系数对应关系

前因子来自本机固定 v2 源 `main.tex:592–604`：

\[
R^{\rm in}=(r-r_+)^{-2iMr_+(\omega-m\Omega_H)/(r_+-r_-)}
\sum_{j=0}^4 A_j(r-r_+)^j,
\]
\[
R^{\rm up}=\underbrace{\frac{e^{ikr_*}r^{iM\mu^2/k}}{\sqrt{r^2+a^2}}}_{P(r)}
\sum_{j=0}^4 B_jr^{-j},\qquad
k=\sqrt{\omega^2-\mu^2},\quad A_0=B_0=1.
\]

`ΩH=a/(2Mr+)`。正文说使用到 j=4；notebook 提供对应完整系数，但不再包含 R(r) 或积分初值代码。其说明位于 `Boundary Coefficients.nb:22–32`，`χ` 是无量纲自旋 `a/M`，`Λ` 是 massive Klein–Gordon 的**角向本征值**，不是附加 `a²ω²−2amω` 后的 Teukolsky 分离常数。

| 系数 | notebook 赋值起始行 |
|---|---:|
| A1 | 43 |
| A2 | 132 |
| A3 | 585 |
| A4 | 1134 |
| B1 | 2290 |
| B2 | 2344 |
| B3 | 2495 |
| B4 | 2873 |

所有八个 Input cell 已做不求值的 Box 语法转换，保留括号，导出 `coefficients.wl` 与带行号的 `coefficients_inputform.json`。仅涉及 RowBox、FractionBox、SuperscriptBox、SqrtBox；没有运行 notebook 内的求解器，因为它不含求解器。

特别是作者 B1 可直接化简为

\[
B_1=-\frac{M\mu^2}{2k^2}
+\frac{i}{2k}\left[\Lambda-4M^2\mu^2+a^2k^2+\frac{M^2\mu^4}{k^2}\right].
\]

这个表达式的 `k^{-3}` 行为清楚显示：临近传播阈值，固定阶数的 `1/r` 展开不一定能在惯常外边界半径内使用。`j=4` 本身不是收敛保证；还必须检查各项及外边界位置。

## 原系数与独立递推的交叉核对

另一个模块从精确 `P'/P` 和径向 ODE 独立递推四阶系数。本次审计单独读取作者 notebook 表达式，用受限算术解析器（不调用 Python eval）与 mpmath **80 位**运算求值。在现有两组 scalar22 参数、传播与束缚两侧，作者 B1..B4 与独立递推结果的最大相对差为 **1.4333×10^-45**，受对照输出的 45 位字符串截断限制。

这是**原作者表达式的独立数值核对**，不是原 Mathematica 程序执行的声明，也不是同一递推自我比较。它确认当前 `Li_exact_prefactor_series4` 控制确实用了公开的同一组 B1..B4；不能再把它与其他前因子下的零阶/六阶边界混称。

| 本地轨道 | ω | 作者 B1 | `|B1|/4000` |
|---|---:|---|---:|
| 41.1 | 0.30007587885635806 | −988.292654128 + 13599.407993700 i | 3.408817795 |
| 42.1 | 0.2999423381255306 | −18124.524456460 | 4.531131114 |

这里代入的是本地同步背景 `a=.8771530275949366, μ=.3, ell=m=2`，并非声称作者计算使用了这些精确参数。首修正项在 `R=4000` 已超过 1，这直接提示该位置不足以凭低阶渐近假设保证边界正确。它不单独给出总场误差，更不能识别作者实际 R。

## 没有公开的运行设置

对唯一 notebook 全文搜索 `NDSolve, NIntegrate, Rout, Rmax, Rmin, RUp, WorkingPrecision, PrecisionGoal, AccuracyGoal, 4000, 1000` 均为零次匹配。更强的结构性证据是全部内容仅一个说明 Cell 加八个系数赋值 Cell，没有径向积分或参数实例化。

所以证据层级是：

- **已证实**：论文前因子、j≤4、A0=B0=1、八个公开系数及其参数含义。
- **已独立核对**：B1..B4 在两组实际近阈值参数下与形式递推一致。
- **仍未知**：原图 Rout、近视界 offset、实际采样/积分容差、是否按轨道动态增大边界。
- **不能推出**：作者原图的差异必由 Rout4000、转折点以内初始化或某个固定截止值导致。

## 产物

- `outputs/paper_original_reference/dongjun826_EMRI_in_scalar_clouds/`：原仓库固定提交。
- `outputs/paper_original_reference/li_boundary_audit_20260917/coefficients.wl`：八个完整作者表达式的可读 Wolfram 形式。
- 同目录 `coefficients_inputform.json`：原行号与表达式。
- 同目录 `verify_author_coefficients.py`：独立、无任意代码求值的 80 位数值核对。
- `docs/field_alignment_20260917/li_author_boundary_coefficients_check.json`：八组核对值、参数和 SHA。
- 本报告配套 `_manifest.json`：来源固定提交、范围与文件哈希。
