# Kerr 点粒子 Teukolsky 符号验证

- 来源 TASK：T006 / W02.2
- 状态：owner-reviewed executable evidence；科学结论和边界见 [findings](../../findings/kerr_point_particle_evolution/README.md)
- Theory targets：[DERIVATION.md](../../findings/kerr_point_particle_evolution/theory/DERIVATION.md)
- Result JSON：[data results](../../data/kerr_point_particle_verification/README.md)
- 环境：Wolfram 14.2.1，MacOSX-ARM64；任务执行通过 Mathematica MCP
- 依赖：Wolfram Language built-ins；未使用 PowerExpand；各文件 assumptions 以返回值和 JSON 为准

## 目录

checks/ 完整保留迁移前的 25 个 WL artifacts：

- A1_A2_operator.wl、A2_four_source_blocks.wl；
- A3/A4 compatibility entries；
- a3/：geodesic、两条 stress-tensor routes、tetrad projections 与 summary；
- a4/：fixed-m theta、direct-y source 与 summary；
- a5/：common jets、四块 source、core 与 summary；
- a6/：first-order、characteristics、axis/projection targeted checks 与 summary。

Compatibility entries 和 helpers 均是已审阅执行链的一部分，不因存在 staged summary 而删除。

## 重放入口

在仓库根目录通过 Mathematica Kernel/MCP 执行：

~~~wl
Get[FileNameJoin[{
  Directory[],
  "code", "kerr_point_particle_verification",
  "checks", "a6", "a6_targeted_summary.wl"
}]]
~~~

A3--A6 summaries 从 data/kerr_point_particle_verification/results/ 读取已保存的 stage JSON，只做轻量汇总，不重新执行重代数。需要复算数学证据时，应直接执行对应 stage WL，再由其 summary 汇总；transport failure 与 mathematical residual 必须分开记录。

## 证据强度

- A2--A4 包含完整或 staged endpoint comparisons；
- A5 是 common-jets 与 generic jet-interface blocks 的模块化闭包；
- A6 只包含 owner 要求的 targeted checks，summary 明确 FullEndpointClaimMade=false；
- A1 的 signature/Weyl、sourced-unit 和 independent PDF transcription 边界不由代码位置改变。

这些 checks 是后续 Python contract tests 的 oracle，不是数值演化 runtime dependency，也不应自动生成或修改生产公式。
