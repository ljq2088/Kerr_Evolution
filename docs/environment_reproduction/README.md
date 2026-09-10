# Kerr 环境扰动论文复现

依据：Dyson et al., [arXiv:2501.09806v1](https://arxiv.org/abs/2501.09806v1)，正文及补充材料。
目标是有质量复标量云被圆赤道轨道小天体的引力度规扰动驱动的响应。

## 已有 s=-2 演化核查

本分支从 `2d690e2` 创建。`code/kerr_point_particle_solver` 已有点粒子源
s=-2 Teukolsky 演化，采用七点有限差分和经典 RK4；详见
[`../upstream_alignment.md`](../upstream_alignment.md)。
已存 smoke 数据在 `../sminus2_smoke/`：1301 个样本，T=0.25 至 325.25。
这些是正则四标架中的 peeling psi4 模态；不是完整 Lorenz 规范 h_mu_nu。
上游 Weyl 符号/归一化 bridge 仍 conditional-open，smoke 不代替连续极限收敛。

## 本次完成的实际计算

`src/environment_cloud.py` 在 M=1 下独立求 alpha=.3 的 |211> 阈值云。
使用球谐矩阵对角化角算子，以双向径向 shooting 求解束缚态谱条件，
并同时满足 omega_c=Omega_H。这不是论文采用的 Leaver 实现。

径向算子为

    (Delta R')' + [K²/Delta - mu²r² - a²omega² + 2am omega - Lambda] R = 0.

阈值处采用视界正则 Frobenius 分支；远方采用
R'/R=-k+beta/r，k=sqrt(mu²-omega²)，beta=(2omega²-mu²)/k-1。
该远方条件是渐近截断，程序扩大远区边界、缩小视界偏移并收紧积分容差核对结果。
当前根区间明确限定 alpha=.3 的 |211> 分支，不宣称通用准束缚态求解器。

结果记录在 `cloud_211.json`：

- a/M = 0.87715302759493；M omega_c = 0.29629324847975。
- m=2 通道开启半径 = 41.6608484179 M，与论文 41.66 M 的报告精度一致。
- 两次边界/容差设置所得半径相差约 2.8e-10 M。

通量接口采用统一 C=1、角函数球面积分模平方为1的约定，剥去 epsilon²q²。
它分别返回荷通量、普通波能流和扣除云能量后的轨道有效通量。
这里输入振幅是接口参数；尚未从完整 Lorenz 源算出物理响应振幅。
有质量场的定频大半径通量也不能等同于无质量场的 scri 时域波形。

执行：在仓库根目录运行

```bash
.venv/bin/python src/environment_cloud.py
.venv/bin/python -m pytest tests/test_environment_cloud.py -q
```

测试使用独立 SciPy 扁球谱、论文报告的传播阈值、边界截断变化、
荷与能量记账、视界面积定律检查。阈值一致不证明源归一化或所有论文曲线已复现。

## 后续完整复现步骤

1. 获得或实现圆轨道完整 Lorenz 度规重构，包含质量、自旋与规范 completion；
   核对源单位、规范条件及低多极。已有 s=-2 曲率场不能直接替代此步骤。
2. 输出并按云质量归一化背景径向/角向本征函数，计算协变 Hessian，
   检查 Box phi0 = mu² phi0。
3. 构造 h^{mu nu} nabla_mu nabla_nu phi0，按
   omega_m=omega_c+(m-m_b)Omega_p 将 Sigma 乘源投影到扁球谐。
4. 解有源 massive 径向方程，验证 Delta W 恒定与 In/Up 两端条件，获得物理振幅。
5. 重建尾迹并做模态、径向与源截断收敛；比较论文通量曲线、阈值与共振。
   在比较绝对通量前核对 v1 的场/荷流归一化和图4 epsilon 幂次差异。

当前交付是复现的第一项数值基准；完整 Lorenz 源、受迫尾迹及通量曲线尚未完成。
