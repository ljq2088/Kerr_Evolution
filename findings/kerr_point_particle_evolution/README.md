# Kerr 点粒子 Teukolsky 演化

- 来源 TASK：T006--T007 / W02
- 审阅日期：2026-09-04
- 状态：owner accepted；连续理论、符号证据和数值 solver 已完成
- 主推导：[DERIVATION.md](theory/DERIVATION.md)
- 约定：[CONVENTIONS.md](theory/CONVENTIONS.md)
- Mathematica 代码：[verification checks](../../code/kerr_point_particle_verification/README.md)
- 机器结果：[verification results](../../data/kerr_point_particle_verification/README.md)
- 状态审计：[audit report](evidence/theory_status_audit.md)
- 数值结果：[numerical evolution](numerical_evolution.md)
- 正式代码：[solver](../../code/kerr_point_particle_solver/README.md)、[公共工具](../../code/kerr_waveform_tools/README.md)
- 正式数据：[SXS0305 medium bundle](../../data/kerr_point_particle_evolution/sxs0305_medium/manifest.json)

## 理论链

~~~text
Sasaki/Kinnersley sourced equation
  -> Ripley field、tetrad与hyperboloidal coordinates
  -> 完整hyperboloidal four-block source
  -> constant-T point-particle stress tensor与三个Kinnersley projections
  -> fixed-m与y=-cos(theta)
  -> exact distributional equation
  -> Gaussian / compactified-Delta extension与解析time jets
  -> LL、LJ、JJ、JL四个smooth source blocks
  -> A--F二阶PDE与P/Q一阶系统
  -> characteristics、axis regularity与SCRI+ projection contract
~~~

正文保留推导主线，theory/appendix/ 保留完整 coefficients、Jacobians、source nesting、time jets、四块数值源和 A6 数值接口。迁移没有用摘要替代这些细节。

## 证据状态

| Gate | 正式证据 | 可支持结论 | 开放边界 |
|---|---|---|---|
| A1 | A1/A2 operator check | 输入的 boost/spin relation 下，field factor 与12-jet operator algebra一致 | signature/Weyl bridge、sourced units、independent Ripley Eq. (22) transcription |
| A2 | complete BL-to-hyperboloidal source check | four-block source 与总 endpoint residual 为零 | primary-PDF transcription 与 A1 |
| A3 | geodesic、constant-T、BL pullback、tetrad stages | 三个 point-particle Kinnersley components 的 local endpoint closure | global monotonicity、unique slice intersection、turning-point differentiability、A1 |
| A4 | staged unprojected-to-fixed-m/y checks | exact distributional field/four-block endpoint residual 为零 | 数值 axis/Gaussian/discretization 与上游 assumptions |
| A5 | common jets 与 four-block interface checks | interior time-algebra interfaces 模块化闭合 | 不是单一完全展开生产源；finite-domain、turning point、数值实现 |
| A6 | first-order、characteristics、axis/projection targeted checks | 三个窄接口 checks 通过 | FullEndpointClaimMade=false；production boundary rows、stability、convergence、transition |

## Owner 审阅边界

Owner 接受当前连续方程、推导链和已记录的符号检查，但明确不认为 horizon、SCRI+ 与 axis 附近的 production finite-difference rows 已被充分理解或验证。A6 的 indicial derivation 和 toy factor-aware basis 只说明局部 regularization 有一致的数学构造，不替代 W02.3 的 boundary-closure note、实际 matrices、manufactured tests 和稳定性证据。

## 文件导航

- [Ripley operator](theory/appendix/RIPLEY_OPERATOR.md)
- [Sasaki source](theory/appendix/SASAKI_SOURCE.md)
- [Point-particle projections](theory/appendix/POINT_PARTICLE_PROJECTIONS.md)
- [Fixed-m 与 y equation](theory/appendix/W02_1_4_FIXED_M_Y_DERIVATION.md)
- [Gaussian extension 与 time jets](theory/appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md)
- [四块数值源公式](theory/appendix/a5_four_block_source_terms.md)
- [A6 连续演化与数值接口](theory/appendix/a6_numerical_interface.md)

## 下游接口

W02 已完成。历史计划继续位于 [workspace](../../workspace/plan/W02_point_particle_evolution/PLAN.md)；正式 solver 从这里读取理论和约定，不建立第二份可编辑公式。Mathematica residual 是可重放证据，不是 runtime dependency。数值结果和下游 W03--W06 接口见 [numerical evolution](numerical_evolution.md)。
