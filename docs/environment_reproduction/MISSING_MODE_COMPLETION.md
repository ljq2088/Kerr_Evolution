# 按清单补齐受迫标量通道

实现：src/report_missing_environment_modes.py；批处理接口：src/report_environment_mode_batch.py。

输入是 report_flux_coverage.py 的清单。模块读取清单指定的响应文件，核对模态、背景参数、分辨率、Green 边界、完成状态和通量记录。缺失文件不会因为旧清单写着完成而被忽略。目标取 field、infinity、horizon 的去重并集，不再从场的 ell>=2 范围推断所有通量已齐。

非静态通道按 |m_g| 和 Green 边界设置分组。只有远边界半径和构造方法一致的正负 m_g 才合并。批处理现支持省略 --scalar-ells，仅指定 --conjugate-ells，复用正 m_g 度规缓存，求解负 m_g 的目标通道；不再需要附带一个多余的正分支通道。原有正分支及混合批次接口保持兼容。

静态 m_g=0 被单独列为 unresolved_static_modes。执行模式遇到静态缺项会在开始计算前拒绝执行，要求独立匹配数据。不会将静态场误当作零，也不会凭空指定补全系数。

默认只生成计划。--execute 依次运行批次并记录返回状态；失败时保留已完成通道，停止后续批次。结束状态明确要求刷新覆盖清单，不将执行成功当成收敛证明。重跑前应重新生成覆盖清单。

```bash
.venv/bin/python src/report_missing_environment_modes.py <coverage.json> --output <plan.json>
.venv/bin/python src/report_missing_environment_modes.py <coverage.json> --output <plan.json> --workers 2 --execute
```

## 实际验证

41.6M 和 41.8M 清单各缺 scalar(ell,m)=(1,-1)。规划器分别只请求 m_g=2 的共轭 ell=1，Green 外边界4000M。两批各读取96个缓存径向点，重构点数均为0；输出均完成96点源积分。批次为 scalar_batch_mg2_9dbfa94712ae.json、scalar_batch_mg2_9fdf1cfdaff4.json。

刷新后两轨道 field=88/88、horizon=85/85；infinity分别36/36、30/30。对刷新后的清单重新规划，缺项均为空。figure5_complete_inventory_audit.json 记录每轨道90个唯一响应文件的指纹和有限值核验。有限总通量为：

| rp/M | 无穷远轨道有效能流 | 视界轨道有效能流 |
|---|---:|---:|
| 41.6 | 3.3379805913895243e-5 | -3.8135563126411016e-6 |
| 41.8 | 1.897384758172336e-5 | -3.7925652627451056e-6 |

单位保持既有的每 q^2(Mc/M) 归一化，没有拟合幅值。

10项测试通过：原有批次失败保留结果、共轭单独运行、空请求拒绝、低阶漏项、不同边界拆分、重复模态去重、静态缺项、清单缺失文件、错误模态和分辨率拒绝。上述检查不替代度规源、径向/角向截断、边界收敛或论文归一化验证。
