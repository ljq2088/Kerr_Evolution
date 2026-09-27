# 径向解缓存修复（2026-09-25）

用户在旁支对话要求修复性能问题，并要求更大的缓存。最终仅把 src/lorenz_kappa.py 的 scalar_resolvent LRU 上限从96改为512；128是准备阶段，未启动128版计算。

当前 L=20、|mg|=1 的完整扫描需要100组辅助径向解。原96项LRU在逐角度遍历时循环淘汰。新测试覆盖两个半径、每半径40个角度：100次未命中、7900次命中。3组真实径向解与无缓存重新积分逐元素一致；既有kappa方程和径向域检查也通过，共21项测试。

改动前后除缓存装饰器及注释外AST完全一致。已保存的87个半径点按明确的缓存策略兼容迁移到新缓存目录，逐文件检查张量数组相同，并保存旧元数据、原文件散列及迁移记录。没有把旧数据冒称为重新计算。旧任务的4个未完成半径可能需要重算。

旧进程206941已停止。新进程269447已启动，保持4个工作进程，物理参数、精度及所有截断不变。

当前执行记录：
`docs/li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512/execution.json`

兼容迁移记录：
`docs/li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512/cache_policy_migration.json`

原目录 execution.json 已标记 superseded_by_cache512_restart，并指向新目录。全局 implementation_manifest.json 的 active_run 已更新。以后查询应读取新目录。

本次未测量整体加速比，不宣称40倍。128已足以容纳当前100项；512为更高角阶提供余量。进程仅按需保存条目，不预分配512份解。后续|mg|>=2时所需条目本来较少，不能按第一组的加速情况外推所有模式。

相关测试：tests/test_lorenz_kappa_cache.py、tests/test_lorenz_kappa.py、tests/test_lorenz_kappa_domain.py。
