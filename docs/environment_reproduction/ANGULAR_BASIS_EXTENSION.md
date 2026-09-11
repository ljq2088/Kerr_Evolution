# 高阶度规加密中的角向基底越界修复

rp=10M、Lg=24、nt=24任务以退出码1结束，异常来自environment_cloud.angular_eigenvalue固定20维矩阵：ell=21,m=1时索引20越界。没有完成任何径向度规采样或标量通道。失败批次单独保存在scalar_batch_mg1_d7551b4b26ac_failed_fixed_basis.json，含原日志SHA256；原日志保持不变。

角向矩阵现随目标ell增长，目标模态以上保留至少11个基函数，并校验量子数与最小基底大小。两端cos(theta)乘积所需的额外一维继续保留。size参数表示最小矩阵维度。

测试覆盖ell=18,20,21,24,32，m=0,1,5，c=.1,1,3的正负c²：与SciPy独立的pro_cv/obl_cv及80维矩阵比较。云、角谱和批次失败测试共6项通过，径向与协变源测试9项通过。另对rp=3.5,10,20,30、所有1<=m<=ell<=18和质量平方偏移0,+-5e-5比较旧新本征值，最大绝对变化1.3074e-12，见angular_basis_extension_audit.json。此检查不是最终度规或通量误差界，特别是kappa使用本征值相关的参数差分。

批次度规预计算异常现在也会写入明确失败状态，避免已退出任务仍显示batch_in_progress。增加了异常路径测试。源码变化产生新的度规缓存哈希，不冒充旧缓存；原有已启动队列继续运行，后续批次使用修复后的代码。

在相同Lg=24、nt=24、nr=8参数下恢复两进程任务，日志outputs/rp10_L24_nt24_mg1_extended_basis.log。该任务同时改变度规与角积分阶数，不是单独的Lg收敛实验。完成前不报告加密通量结果。

独立高阶径向预检查随后发现environment_source.angular_mode也使用固定20维基底。为避免已知确定性越界，核对PID和进程组后终止了首次恢复任务，保留单独的stopped_eigenfunction_basis记录。这不是因观察超时而重启。现在本征值与角函数共用angular_basis_size；新增高阶角函数归一化、导数加密和径向本征值一致性测试，相关测试共16项通过。ell=21..24各取质量平方0,+-step,+-step/2的20个实际径向Green系统通过初始化，且ell=24,m=1在r=5,theta=1.2的完整重构度规数值有限，见high_ell_resolvent_preflight.json。这些是初始化与有限性核验，并非高阶度规收敛证明。再次恢复的日志为outputs/rp10_L24_nt24_mg1_shared_basis.log。
