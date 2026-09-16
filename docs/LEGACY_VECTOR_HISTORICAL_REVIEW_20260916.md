# 公开 spin-1 下界错误与论文图 2 的历史归因审查

审查日期：2026-09-16。范围为用户提供的 `Kerr_environment_horizon_flux_diagnosis.md`、固定 GitHub 提交及原论文 arXiv v1 TeX；本轮未修改物理实现、未做新的数值消融，也未联系作者。输入文件的 SHA256 及全部取证文件见随附 manifest。

**结论：公开代码的错误和修正内容可以确认；它与 2025 年论文图 2 的实际运行版本之间的对应关系尚不能确认。** 用户诊断第 14 节把“论文使用模块化 Mathematica 包”进一步写成“图 2 使用了这个有缺陷的公开生成路径”，超出了目前原始证据。另一个必须纠正的论点是：这个门槛不含 Kerr 自旋参数，因此它也会删除同一路径下 Schwarzschild 的相应 spin-1 模态。不能仅凭 Schwarzschild 是非旋转背景就排除它作为反证测试。

## 1. 提交确实包含什么

[官方提交记录](https://github.com/ConorDyson/KerrLorenzMSF/commit/b784aa25f7d7a226389e2a5fbdac1a883fd411ff)及其 [REST 元数据](https://api.github.com/repos/ConorDyson/KerrLorenzMSF/commits/b784aa25f7d7a226389e2a5fbdac1a883fd411ff)给出：

- 提交：`b784aa25f7d7a226389e2a5fbdac1a883fd411ff`。
- 作者及提交时间均为 **2026-06-02 22:51:59 UTC**。
- 提交说明：`additional ll>lmax2 bug found`。
- 父提交：`d7d7f829bfb57b2e6263da3a0ce1207f499c1beb`，时间 2026-05-21 14:19:50 UTC。
- 提交统计为 796801 行增加、0 行删除。

因此，提交说明本身不能作为原路径被就地修复的证据。完整递归 tree 均未截断；父树与新树的逐路径 blob 比较结果为 **1888 个新增文件、0 个删除文件、0 个既有文件内容修改**。这是一次批量加入新目录的提交。

但是，以下**不同路径之间**的内容差异确实直接证实了修正：

| 版本与文件 | spin-1 下界 | Git blob SHA |
|---|---|---|
| 父提交 `GenerationCodes/Numerics-WorldTube/h1-radial-gen/MetricReconstructRadiative.wl` | 错误的 `lmins2` | `deec703c45fb5f7e6e5a8630c33cb7eb4caaa2f6` |
| b784aa 新增的 `Ret-Numerics/h1-radial-gen/MetricReconstructRadiative.wl` | 正确的 `lmins1` | `4be336657bda0a44dbc44b77b2b397725d10072d` |

[父提交原文件](https://raw.githubusercontent.com/ConorDyson/KerrLorenzMSF/d7d7f829bfb57b2e6263da3a0ce1207f499c1beb/GenerationCodes/Numerics-WorldTube/h1-radial-gen/MetricReconstructRadiative.wl)与[新增修正版](https://raw.githubusercontent.com/ConorDyson/KerrLorenzMSF/b784aa25f7d7a226389e2a5fbdac1a883fd411ff/Ret-Numerics/h1-radial-gen/MetricReconstructRadiative.wl)的完整差异只有两个门槛替换及一个空行删除：

1. 旧文件第 2335 行／新文件第 2334 行：生成 `MMm1vec` 中 `s=-1` 的 `MSTMode`，`If[ll>=lmins2,...]` 改成 `If[ll>=lmins1,...]`。
2. 旧文件第 2431 行／新文件第 2430 行：构造 `MMm1gridL,MMm1gridR`，同一门槛改成 `lmins1`。

完整小差异文件：[LEGACY_VECTOR_CROSS_PATH_20260916.diff](/home/ljq/code/kerr-hyperboloidal/docs/LEGACY_VECTOR_CROSS_PATH_20260916.diff)。

**同一提交中仍保留有错误版本。** 旧 WorldTube 路径未修改；新加入的 `Ret-Numerics/BIGRUN-h1-radial-gen-session/MetricReconstructRadiative.wl` 也使用旧文件同一个错误 blob。因此准确表述是“该提交加入了修正版副本”，而非“该提交修复了仓库中所有生成路径”。提交标题里的 `lmax2` 也不是对真实差异中变量 `lmins2` 的精确描述。

## 2. 公开历史不能锁定论文运行版本

固定到 b784aa 的路径历史显示，旧 WorldTube 文件在 2026-05-04 的 `c34221d1f8acc2924b7a4ecf04ecf30fd1444680` 和 2026-05-21 的父提交中有记录；新 Ret 文件在 2026-06-02 首次出现在该路径下。已保存的 GitHub 仓库元数据显示仓库创建于 2026-01-15。**这些是公开仓库记录时间，并不是文件首次编写或数据实际计算的时间。** 它们既不能证明代码在 2025 年已存在并被使用，也不能证明它当时不存在。

论文 [arXiv v1](https://arxiv.org/abs/2501.09806)发表于 2025-01-16；本次逐行审查的 [原始 TeX](https://arxiv.org/src/2501.09806v1)中：

| 原始位置 | 能支持的事实 | 不能支持的推断 |
|---|---|---|
| `main_PRL.tex:681`，补充材料 S.5 Metric data | 作者将早期 Mathematica notebook 改成适合大参数空间的模块化包，输出 spin-weighted spherical 数据，再组合到二维网格 | 未提供该包的公开路径、代码散列或与 2026 仓库的逐文件映射 |
| `main_PRL.tex:343–344` | 图 2 使用 `Figures/Flux_Inf_Hor_TwoPanels_2.pdf`，展示总标量通量，视界通量求和至标量 ℓ=5 | 图文件不是运行日志，不能反推当时是否经过错误门槛 |
| `main_PRL.tex:362` | 论文报告视界通量主要由标量 (ℓ,m)=(0,0) 贡献，且有已说明的高半径例外 | 标量单极主导并不单独证明度规 ℓ=1 vector 被删，或证明其影响的符号、大小 |
| `main_PRL.tex:685` | Schwarzschild 极限下度规数据经过独立 Lorenz 数据检验 | 不能把这段改写为“Schwarzschild 生产一定走另一套绕开 bug 的代码” |

原 TeX 没有 GitHub 地址、commit、图 2 生成脚本散列或相应原始中间数组。此前公开材料检索也尚未获得与图 2 一一对应的完整环境程序与原始通量表。因此诊断第 14 节应改为：

> 论文 S.5 确实说明使用模块化 Mathematica 包。2026 年公开代码中存在可识别的门槛错误及修正版，但目前没有确认它们与 2025 年图 2 运行版本的对应关系。

公开 HDF5 样本更接近某个删项版本，若数值比较可靠，可支持“该样本可能受相同问题影响”；它不自动成为图 2 原始流水线的版本证据。尤其应分别记录样本参数、生成来源和图 2 参数，不能把一次单模态样本的身份扩展为整套论文结果的身份。

## 3. 为什么 Schwarzschild 对这个假说仍有检验力

旧 WorldTube 文件第 2298–2300 行定义：

```wolfram
lmins0 = Abs[mm];
lmins1 = Max[1, Abs[mm]];
lmins2 = Max[2, Abs[mm]];
```

对 `|m|=1, ell=1`，正确条件是 `1>=1`，错误条件是 `1>=2`。两者的区别完全不依赖 `a`。同一函数第 2288–2296 行还明确处理 `a0=0`，随后进入上述门槛；并没有一个 Schwarzschild 特例绕过两个错误的判断。

因此，在**同一个错误生成路径**中，`a=0` 也排除了 `ell=1, |m|=1` 的 spin-1 径向构造。对非零 `m` 的圆轨道，Schwarzschild 轨道频率仍非零，不能自动归入静态 `m=0` 的另一条路径。

这只严格证明“模态被删除”，并不未经计算就证明其最终通量改变量一定非零：特定几何、源项或观测量有可能抵消。若要主张该删除只改变 Kerr 结果，必须给出相应 Schwarzschild 向量源为零或贡献相消的推导，或给出同路径消融结果。**诊断第 12 节仅以非旋转极限为理由排除这个检验，在逻辑上不足。**

而且论文 S.5 的验证叙述提供了一项具体约束：作者称 Schwarzschild 度规的全部模式和半径与独立 Lorenz 数据达到机器精度，仅 `(ell,m)=(1,0)` 因 completion 选择有差别，并在加修正后消除。这里的例外不是门槛假说涉及的 `(1,+1)` 或 `(1,-1)`。这一点与“论文实际用了直接丢失该 vector 的同一路径”之间存在需要解释的张力，不能把两种低模态问题混为一谈。

论文关于 Schwarzschild 云并非严格定常、手动令云频率虚部为零的讨论（TeX:382），以及改用 ingoing Eddington–Finkelstein 坐标重算的验证（TeX:697），涉及背景云、导数和坐标正则性；它们不是这个自旋无关门槛在 Schwarzschild 中不起作用的证据。

## 4. 证据等级与可被检验的结论

| 判断 | 当前等级 |
|---|---|
| 公开模块化代码存在错误的 spin-1 最小 ℓ 门槛 | **已证实：原代码** |
| b784aa 加入的特定 Ret 副本把两处 `lmins2` 改成 `lmins1` | **已证实：固定提交跨路径差异** |
| 同提交所有生成路径都已修复 | **不成立：错误副本仍存在** |
| 门槛错误只影响 `a≠0`，所以 Schwarzschild 不能检验 | **不成立：门槛和分支代码不支持** |
| 论文使用过某种模块化 Mathematica 度规数据生成流程 | **已证实：论文 S.5** |
| 该流程就是公开仓库里的这个错误副本 | **尚未证实：缺历史对应关系** |
| 这个错误解释当前 17–33% 的 Kerr 视界通量差距 | **需数值检验：不是代码检查本身的结论** |
| 删项后拟合图 2 就证明作者当年确实用了它 | **不足：定量支持候选解释，不等于历史版本取证** |

在线性响应中，如果删除 vector 导致振幅由 `Z_full` 变为 `Z_full-δZ_vector`，则通量变化含干涉项：

```text
F_legacy - F_full = C * ( -2 Re[conj(Z_full) δZ_vector] + |δZ_vector|^2 ).
```

因此，指出缺失通道和角动量选择规则，尚不能确定通量偏差的符号或百分比。后续消融应明确使用相同背景、轨道、云、投影、径向 Green 函数和归一化，仅切换该 vector 通道，并在 Kerr 与 Schwarzschild 的多个参数点比较。若变化幅度或方向不对，它便不能解释相应偏差；即使定量一致，也应保留“候选历史解释”的限定，直到有版本或运行数据证据。

本审查没有据此断言论文结果错误。当前最强可复核结论是：**存在一个真实、可隔离测试的公开代码问题，但原诊断对其论文历史归因和 Kerr 专属性的表述过强。**

## 5. 取证记录

- 小 manifest：[LEGACY_VECTOR_HISTORICAL_REVIEW_20260916_manifest.json](/home/ljq/code/kerr-hyperboloidal/docs/LEGACY_VECTOR_HISTORICAL_REVIEW_20260916_manifest.json)。记录官方 URL、固定 SHA、取证文件 SHA256、原 TeX 行号与本轮范围。
- 原始 API、完整源文件和输入副本保存在 `/home/ljq/code/kerr-hyperboloidal/outputs/paper_metric_reference/legacy_history_20260916/`。GitHub commit 的默认文件列表存在分页限制，本结论采用未截断的递归 tree 完成全提交核验。
- 本轮新增审查文档与小 diff；没有改动生产代码或历史输入文档。


## 6. 追加：两个门槛的因果范围与删除实验的等价条件

本节只读检查同一固定原文件的数据流，未重跑 Mathematica。结论是：**在已核对的两个副本中，这两处错误不会通过 spin-1/κ 的联立匹配反向改变 κ、trace 或更高 ℓ 的系数。** 它们发生在径向模式求值和构造数组阶段，匹配系统仍完整保留最低 spin-1 模态。

以下行号均指保存的父提交 WorldTube 文件；修正版除了前述门槛和空行外内容相同。

| 环节 | 原文件行号 | 数据流检查 |
|---|---|---|
| 联立匹配函数参数 | 594–606 | `CalcJumpsκs1` 接收轨道、spin-2/trace 跳跃、源投影、角向矩阵和本征值；不接收 `MMm1vec` 或 `MMm1grid` |
| 匹配内部最低模态 | 679–710 | 内部重新定义 `lmins1`，构造 `Pm1/Pp1` 和 spin-1 跳跃时已使用正确下界 |
| 匹配未知量 | 979 | `Pm1j0/Pm1j1` 均从 `lmins1` 开始，κ 从 `lmins0` 开始 |
| `|m|=1` 额外闭合条件 | 1009–1015 | 保留两条 `ell=1` 条件，没有因为 MST 数组为零而删除未知量或重解低维系统 |
| 联立求解 | 1018–1022；主函数 2390–2403 | 一次性求出完整 `s1κjumps`，后续不重写这些跳跃 |
| 错误径向门槛 | 2335；2429–2437 | 前者跳过 `ell=1` 的 MST 对象，后者直接返回两侧的 `{0,0,0}` 径向数组 |
| 生成正自旋函数 | 243–249；2446–2447 | `Pp1fromPm1` 对输入径向函数及导数是线性的；零输入产生零的另一手性 |
| κ/trace 径向解 | 75–87；2450–2458 | 只使用 `MMhκvec`、trace 跳跃和已经完整求得的 κ 跳跃，不读取被置零的 spin-1 数组 |
| 角向重构及总和 | 1895–1924；2461–2476 | 分别计算三个扇区，最后直接 `sall=s2+s1+s0`；没有向匹配或 spin-0 分支反馈 |

因此，对于这里的 `|m|=1` 非静态模式，若固定同一个生成器的其余数据和匹配解，两个门槛在数学上的效应是

```text
h_bad = h_correct - h_vector_from_spheroidal_ell1 .
```

这不是删除联立系统中的 spin-1 未知量后重新匹配：那样会改变 κ 等系数，属于不同实验。该生成器把 κ 的 trace 驱动部分和齐次匹配部分都保留了；其角向非对角修正也只依赖 trace、Γ 矩阵和 spin-0 本征值（691–703、1932 起、2117–2135），不读取 `MMm1grid`。代码中没有发现由这两个数组触发的独立 χ 重匹配分支。

还必须区分**被截掉的 spheroidal 径向标签 ℓ=1**与最终输出的**spherical 度规标签 ℓ=1**。Kerr 下 `Bmat` 和角向微分/混合矩阵会把前者的 vector 贡献分配到多个输出球谐分量。等价删除应减去该径向模式生成的完整 vector 张量及其角向后代，而不是把最终输出某一个 spherical ℓ 分量归零。若只删除最终 spherical ℓ=1，不能声称等价。

本地 `report_legacy_vector_ablation.py:20–24` 重构一个 `ell_g=1,m_g=+1` 的完整 spin-1 张量，再共轭得到 `m_g=-1`；其 43、47–54 行在固定云、其他度规扇区、Green 核与通量前因子下线性减去该源贡献。从操作形式看，它针对的正是上述“保留匹配系数后删除 vector”机制。不过，本节的源码审查**没有额外证明本地 `spin1_metric(...full_current=True)` 与原生成器在全部参数、截断、投影及数值细节上完全相同**，也没有将整套原环境程序重新执行。两者之间的数值对应应引用独立的模块对照证据。

对父报告建议使用的范围表述是：

> 已核对原代码，两个错误门槛不反馈修改 κ/trace 匹配；其直接效应是删除最低 spheroidal spin-1 径向模式的完整 vector 贡献。当前同网格删项实验检验了这一具体机制对所选标量通道的影响；它不是对旧生成器、完整环境求解流程或论文历史版本的全面复刻。

单个 `alpha=0.3,rp=20,scalar(ell,m)=(0,0)` 的结果仍不涵盖总通量的其他标量模式、其他参数点或其他潜在旧实现差异。若该实验不支持所需幅度/方向，可反对“这两个门槛通过该被删通道解释该点偏差”的假说；应避免把它扩大为排除所有旧代码问题。
