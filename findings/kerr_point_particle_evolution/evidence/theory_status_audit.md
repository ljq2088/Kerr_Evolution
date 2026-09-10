# W02 理论与符号证据状态审计

> 本文是 T006 结束时的状态快照；其中“当前 workspace”与“待迁移”措辞不表示
> T007 后的仓库状态。T007 的最终入口见[数值演化 finding](../numerical_evolution.md)。

- 范围：W02.1 主推导及全部附录、W02 三份子计划、workspace 导航，以及 A1--A6 当前可见的 Mathematica 结果
- 状态：任务 12 remediation 已同步；owner migration decision 与开放边界处理待定
- 审计方式：静态比对文档声明、可重放代码结构和已保存 JSON；本审计未执行 Mathematica，也未修复任何发现
- 当前结论：任务 12 已同步 A3 附录、A6 状态、验证方法总览、附录 E 时态和 Gate B2 措辞；A5 的验证仍是模块化闭包而非可直接移植的完全展开源。Owner 接受当前连续/符号范围，但 boundary-adjacent finite-difference closure 与 A1 仍开放，不能自动关闭 Gate A。

## 高优先级发现

### 1. A3 附录仍明确否认已经存在的 A3 machine closure【任务 12 已解决】

[附录 C](../theory/appendix/POINT_PARTICLE_PROJECTIONS.md#c1-证据层级与-a3-起点) 顶部仍写“等待新的 A3 Mathematica closure，不宣称 A3 已通过”，末节又把两条验证路线写成将来任务，并在结尾声称 T006 没有执行 A3。这与以下证据直接冲突：

- [A3 endpoint summary](../../../data/kerr_point_particle_verification/results/a3/a3_endpoint_summary.json) 的四项 stage flags 均为 `true`，且 `A3FullEndpointVerified=true`；
- [DERIVATION Section III](../theory/DERIVATION.md#iii-赤道-kerr-点粒子源) 已把 A3 记为 local algebraic closure；
- [workspace README](../README.md#verification-状态) 同样把 A3 记为已验证并继承 A1。

精确位置：`theory/appendix/POINT_PARTICLE_PROJECTIONS.md:6,508-526`；机器证据：`mma/results/a3/a3_endpoint_summary.json:16-35`。任务 12 已把附录状态和 C.7 改为已执行的 local algebraic closure，并保留 global worldline/turning-point 与 A1 边界。

## 中优先级发现

### 2. 审计时 A6 targeted checks 已执行，但导航和计划仍写成“未执行”【本任务后续已解决】

审计落盘期间已经形成：

- [A6 first-order result](../../../data/kerr_point_particle_verification/results/a6/a6_first_order.json)：$A$--$F$、二阶方程恢复、$P$ 定义、$Q$ 恒等式 residual 均为 $0$；
- [A6 characteristics result](../../../data/kerr_point_particle_verification/results/a6/a6_characteristics.json)：径向、角向特征多项式及 SCRI+、视界极限 residual 均为 $0$；
- [A6 axis/projection result](../../../data/kerr_point_particle_verification/results/a6/a6_axis_projection.json)：indicial/Laurent、toy factor-aware basis 和两组 analytic-mode projection 均通过；
- [A6 targeted summary](../../../data/kerr_point_particle_verification/results/a6/a6_targeted_summary.json)：三个窄标志与 `AllTargetedChecksPassed` 均为 `true`，并明确 `FullEndpointClaimMade=false`。

初次审计时，[README](../README.md#verification-状态)、[主计划](../../../workspace/plan/W02_point_particle_evolution/PLAN.md)、[MMA 计划](../../../workspace/plan/W02_point_particle_evolution/mathematica_verification/PLAN.md)、[A6 附录](../theory/appendix/a6_numerical_interface.md) 和 DERIVATION 的 A6 状态仍分别写“尚未推导或执行”“A6 未执行”或“等待最小 checks”。任务 11 后续已统一改为“theory contract drafted + targeted checks passed + owner review open”，并明确没有 full-endpoint claim；README 的理论范围也已同步为 W02.1.1--1.6。

最终状态位置：`README.md:5-6,29`，`theory/DERIVATION.md:4,1258,1265`，`theory/appendix/a6_numerical_interface.md:7,567`，`plan/W02_point_particle_evolution/PLAN.md:6`，`plan/W02_point_particle_evolution/theory/PLAN.md:7`，`plan/W02_point_particle_evolution/mathematica_verification/PLAN.md:7,30`。

### 3. 审计时 `CharacteristicsChecked` 尚未覆盖 angular characteristic【本任务后续已解决】

[A6 验证计划](../../../workspace/plan/W02_point_particle_evolution/mathematica_verification/PLAN.md#a6cprincipal-symbol-与-continuum-boundaries) 要比较

$$
v_R^\pm=\frac{B/2\pm\sqrt{B^2/4-AC}}{A},
\qquad
|v_y|=\sqrt{\frac{1-y^2}{A}}.
$$

初次审计时，[a6_characteristics.wl](../../../code/kerr_point_particle_verification/checks/a6/a6_characteristics.wl) 只把 `(1-y^2)/A` 作为 `AngularSpeedSquared` 输出，未把它纳入 `CharacteristicsChecked`。任务 11 后续已增加

$$
A v_y^2-(1-y^2)=0,
\qquad
v_y^\pm=\pm\sqrt{\frac{1-y^2}{A}},
$$

两个根的 residual 均为 $0$，并已纳入 `CharacteristicsChecked`。更新后的 targeted summary 重放仍给出三个窄标志和 `AllTargetedChecksPassed=true`，同时保留 `FullEndpointClaimMade=false`。

最终证据位置：`mma/checks/a6/a6_characteristics.wl:52-55,74-84`，`mma/results/a6/a6_characteristics.json:4-23`，`mma/results/a6/a6_targeted_summary.json:4-28`。

### 4. A5 的“完整端点”是接口级模块化证明，不是单一完全展开源

A5 common-jets 检查从 Kerr/worldline scalars 验证实际 jets；四块检查则把一般函数 `ff(T,R,y)` 的时间导数替换为抽象 `fJet0,fJet1,fJet2`，再验证 A4 block 与因子化 $\mathfrak H_A$ endpoint 一致。最后 summary 合并两组 flags。

这是一条可成立的模块化闭包：若 common-jets 输出严格满足四块的 jet interface，则四块恒等式适用于实际源。它没有证明失败。但当前证据没有保存“把实际 $A_0,\widehat N_p,\bar M_p,G_R,G_y$ 全部代入四块后得到的单一展开表达式”，也没有生成可直接翻译为 Python 的生产源。因此 `A5FullEndpointVerifiedWithinInteriorScope=True` 应理解为“interior time-algebra interfaces 逐项闭合”，不宜解释成“完整数值实现或完整展开 source 已验证”。

精确位置：`mma/checks/a5/a5_block_core.wl:21-44,76-101`，`mma/results/a5/a5_endpoint_summary.json:10-29`，`theory/appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md:1107-1156`。

### 5. 主计划仍保留已被后续计划替代的 `fixed-width discretization` 验收措辞【任务 12 已解决】

[W02 主计划 Gate B2](../../../workspace/plan/W02_point_particle_evolution/PLAN.md#gates) 仍要求 `fixed-width discretization` 与 Gaussian point-source limit。Theory A6 和 numerical plan 已改为

$$
\sigma_R=4\Delta R,
\qquad
\sigma_y=4\Delta y,
$$

并把默认验收定义为 coupled grid--Gaussian refinement；它明确不是 fixed-source stencil-order 测试。主计划的旧措辞可能重新引入已经放弃的独立 fixed-width gate。

精确位置：`plan/W02_point_particle_evolution/PLAN.md:101-105`；当前接口：`plan/W02_point_particle_evolution/numerical/PLAN.md:94-102,193-203`。任务 12 已把主计划 Gate B2 改为四-cell coupled grid--Gaussian convergence，并明确不单独解释为 fixed-width stencil order。

## 状态矩阵

| Gate | 理论状态 | 当前机器证据 | 最窄可支持结论 | 仍开放 |
|---|---|---|---|---|
| A1 | target/convention candidate | 12-jet operator residual 与输入的 field-factor algebra 为 $0$ | 给定已输入 bridge 与 target transcription 时，算符代数一致 | signature/Weyl bridge、sourced units、Ripley Eq. (22) 独立 PDF transcription |
| A2 | four-block sourced equation 已写全 | BL four blocks 到 hyperboloidal endpoint 的四块及总 residual 为 $0$ | A2 algebraic endpoint verified，继承 A1 | primary-PDF source transcription 与 A1 边界 |
| A3 | local point-particle endpoint 已写全 | geodesic、constant-$T$、BL pullback、tetrad projections 均通过 | 三个 Kinnersley components 的 local endpoint verified | global monotonicity、global unique slice intersection、turning-point differentiability、A1 |
| A4 | fixed-$m$/$y$ exact distributional equation 已写全 | $E_0\to E_m^{(\theta)}\to E_m^{(y)}$、四块和总 residual 为 $0$ | interior exact-distribution endpoint verified | $y=\pm1$ 数值 closure、Gaussian/time jets/discretization、A1--A3 assumptions |
| A5 | interior smooth-source factorization 已写全 | common jets 与 LL/LJ/JJ/JL interfaces 分块通过；total residual 为 $0$ | fixed-width interior time algebra 的模块化 endpoint verified | finite-domain normalization、axis completion、turning branch、coupled point-source convergence、生产代码 |
| A6 | theory-to-code contract draft | first-order、radial/angular characteristics、axis/projection 三组 targeted checks 通过 | 三个窄接口检查通过并已同步；不是 full numerical endpoint | owner contract review、全部 W02.3 工程项 |

## 低优先级一致性与维护发现

- [DERIVATION Section IV](../theory/DERIVATION.md#iv-代数核验与开放边界) 的验证方法总览已在任务 12 更新为 A1--A6 的对象相关 exact-reduction 路线，并明确 A4/A5/A6 的不同证据结构。
- [附录 E](../theory/appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#e11-历史实现的使用与证据强度) 的“未来 A5”时态已在任务 12 改为当前已执行的 A5 checks 与未来 W02.3 tests。
- 已抽查主文、全部附录和三份 W02 计划中的本地 Markdown targets；没有发现缺失的目标文件。现有主要重复是必要的公式摘要与附录展开，而不是第二套互相竞争的权威推导。最明显的过时重复仍是附录 C 的 A3 预验证状态。

## 证据边界详表

### A1--A2

[A1 result](../../../data/kerr_point_particle_verification/results/A1_A2_operator.json) 明确把 `SignatureBridgeAlgebraicallyDecidable` 和 `SourcedUnitBridgeAlgebraicallyDecidable` 记为 `false`，同时给出 12 个 operator jets 全零。因此“算符代数通过”和“A1 物理 convention gate 开放”是相容的；不得由 `OverallSourceFactorResidual=0` 推断单位制已闭合。

[A2 result](../../../data/kerr_point_particle_verification/results/A2_four_source_blocks.json) 记录生成路线未使用 hatted operators、理论 derivative map 或 target block template，四块和总 endpoint residual 均为 $0$。它也明确把 independent primary-PDF transcription 记为 `false`。DERIVATION、README 和计划对 A2 的条件状态与该证据一致。

### A3--A4

A3 summary 没有重新运行重代数，而是汇总四个已保存 stage results；各 stage JSON 均记录 `saved_file_executed=true`，并保留 local-root/global-worldline 边界。因此 summary 的 local endpoint 结论有证据支持。问题只在附录 C 没有同步。

A4 新 staged evidence 明确排除旧 arbitrary-projection template 作为 pass evidence，并从完整周期点粒子 $E_0$ 到 $y$ endpoint 比较。旧 monolithic combined replay 超时，但分组 replay 和新 staged summary 均成功；文档没有把超时误写成非零数学 residual。

### A5

A5 summary 正确保留四项开放边界，并未声称 finite-domain/axis/turning/width 已通过。需要保持的区分是：MMA 证明了 common-jets 模块和 formal-spatial-derivative blocks 的组合接口；空间离散、endpoint rows 及 Python array implementation 尚无机器证据。

### A6

A6 first-order check 支持 $A$--$F$ regrouping、$P/Q$ 连续恒等式和 source 只进入 $P_{,T}$。Characteristics check 支持 radial polynomial、两个 radial boundaries，以及 angular principal polynomial 的两根。Axis/projection check 支持解析 indicial powers、toy rational open-node rows，以及指定 analytic modes 上的连续/离散 Gram matrix。Targeted summary 正确给出 `FullEndpointClaimMade=false`；三个窄标志不覆盖 production rows、quadrature、stability 或 convergence。

## 需要 owner 领域判断的事项

- A1 的 signature/Weyl-scalar bridge、sourced-unit normalization，以及是否接受独立 PDF transcription 的完成方式；这些会影响整体符号或归一化。
- Owner 此前曾指出附录 E “存在公式问题”，随后在任务 12 明确接受当前连续理论与符号范围；当前没有仍待定位的具体 A5 公式编号。迁移时仍须保留 A5 是模块化 interface closure、不是完全展开生产源这一证据边界。
- Post-ISCO transition prescription、$q_{\rm tr}$、flux truncation 与 $(r_0,E_{\rm p},L_{z,\rm p})$ 的准入。
- 是否接受 open-$y$ grid 加局部 factor-aware regularization 作为 baseline，以及 production nodes、axis finite-limit rows 和 quadrature 的数值证据。
- A5 finite-domain Gaussian 处理、turning-point parameterization，以及四-cell coupled grid--Gaussian convergence 是否足以支持生产准入。
- A6 连续 contract 与 targeted checks 完成后是否足以关闭 Gate A。Gate B1/B2 的轨迹 manifest、稳定性、收敛性和 frequency-domain benchmark 不应由 A6 symbolic checks 代替。

## 审计结论

当前 workspace 尚未晋升正式目录。Owner 已接受 W02.1/W02.2 当前连续理论与符号证据范围，但明确不认为 boundary-adjacent finite-difference closure 已完整理解或通过；A1 也仍开放。任务 12 已修复本报告的直接状态冲突，并另写迁移计划。A5 模块化证据边界与全部 W02.3 工程验证继续保留，Full Gate A 仍开放。
