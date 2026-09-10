# Kerr 点粒子 Teukolsky 符号验证结果

- 来源 TASK：T006 / W02.2
- 状态：owner-reviewed machine evidence；开放边界见 [findings](../../findings/kerr_point_particle_evolution/README.md)
- 生成代码：[Mathematica checks](../../code/kerr_point_particle_verification/README.md)
- 环境：Wolfram 14.2.1，MacOSX-ARM64；Mathematica MCP

results/ 完整保留迁移前的 22 个 JSON artifacts，包括 A1/A2 顶层结果以及 A3--A6 staged results 和 summaries。

JSON 中出现的 workspace 路径属于原始执行 provenance，不是当前文件位置。迁移保持数学 residual、assumptions、timing、transport/debug history 和 pass flags 不变；正式导航由本 README 和 findings 提供。

关键 summaries：

- [A3](results/a3/a3_endpoint_summary.json)
- [A4](results/a4/a4_endpoint_summary.json)
- [A5](results/a5/a5_endpoint_summary.json)
- [A6 targeted](results/a6/a6_targeted_summary.json)

A6 的 AllTargetedChecksPassed=true 与 FullEndpointClaimMade=false 必须同时解释。任何 JSON pass 都不替代 owner review、A1 物理边界或 W02.3 的数值稳定性与收敛性。
