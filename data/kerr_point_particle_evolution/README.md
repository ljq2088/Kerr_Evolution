# Kerr 点粒子数值数据

- 来源任务：T007 / W02.3
- 状态：owner accepted
- 数据规范：[数值演化 finding](../../findings/kerr_point_particle_evolution/numerical_evolution.md)
- 生成代码：[Kerr 点粒子 solver](../../code/kerr_point_particle_solver/README.md)

## 正式结果

[SXS:BBH:0305 medium bundle](sxs0305_medium/README.md) 保存 $(2,2)$ 与 $(4,4)$
两套完整结果，包括 complex $\psi_{4,\ell m}$、FFI $H_{\ell m}$、实部/绝对值图、
输入、运行 metadata 和逐文件 SHA256。

## 新计算

正式 CLI 默认先写入 `runs/.staging/<run_id>/`。只有请求的数组、图、metadata、
hash 与 loader 检查全部完成后，才在本目录内原子定稿为 `runs/<run_id>/`。
`status=complete` 与 owner 的 `review_status` 相互独立；失败或中断的运行留在 staging
供恢复或明确清理。首次出现的新背景参数在 `initial_data/` 生成 versioned 初值记录，
相同 key 后续直接复用。
