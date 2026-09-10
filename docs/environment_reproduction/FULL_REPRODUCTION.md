# 完整复现目标与当前证据

目标保持为完整复现 2501.09806v1 正文及补充材料的数值结果。
目前状态：进行中，不能由阈值一致或下列模块测试推断论文已复现。

## 图级验收

| 目标 | 所需证据 | 当前状态 |
|---|---|---|
| 图1，alpha=.3，rp=3.5，赤道与方位切面尾迹，l>=2 | 完整 Lorenz 点粒子源、受迫场、坐标约定、收敛图 | 未完成 |
| 图2，alpha=.2/.3，Kerr 与 Schwarzschild 两端通量 | 正确云归一化，l<=6/5 模和，逐轨道输出与边界收敛 | 未完成 |
| 图3，相对通量误差及旧 Schwarzschild 实现比较 | 同一归一化及规范处理，原始曲线/独立求解器核对 | 未完成 |
| 图4，标量/GW 通量比 | q=1e-6、eta=.1，与作者确认或解释 v1 epsilon 幂次差异 | 未完成 |
| 图5，rp=41.6/41.8 两幅尾迹 | 独立求解两轨道，radiative/bound 过渡，无伪造源 | 阈值位置已核对；场图未完成 |
| 图6，rp=20 的逐模通量 | 模态奇偶与高 l 收敛、主/次模幅度 | 未完成 |
| 图7，rp=20 粒子处 l^-2 收敛 | 点粒子真实源、场连续性、高 l 模和 | 未完成 |
| 共振解释与位置 | 束缚谱和有源响应，普通/轨道有效能流区分 | 未完成 |

## 已实现且验证的模块

1. alpha=.3 的 Kerr |211> 阈值谱：独立双向 shooting。
2. 完整背景 R、R' 与扁球角函数，归一化到 Killing 质量1。
   原文 phi^(1,0) 对应将此场乘 alpha^-3，再由 epsilon=alpha³ sqrt(eta) 得到物理场。
   `background_validation.json`：直接 E-omega Q 约 -2.3e-12；
   Hessian 收缩的最大相对 KG 残差约 3.4e-13。
3. BL Kerr 协变 Hessian 和真实外部 Lorenz h 的收缩/扁球投影接口。
   以有限差分独立核对 Hessian。测试中的 h=g 只是测试张量，不是小天体源。
4. massive 径向 In/Up Green 求解：视界4阶 Frobenius、远方6阶渐近修正。
   正/负频传播和束缚频率的已知解测试；源分段边界放在网格上。
   精确 k=0 被拒绝；近阈值和窄共振仍需专门处理和生产精度验证。
5. Lorenz 重构的 trace h：直接解 Box h=16 pi T。
   已对照公开 2406.12510v3 的 a=.6,ell=m=2,rp=4/6/10/20 振幅表，
   两端复振幅最大相对差约2.4e-9。详见 `lorenz_trace_validation.json`。
   这是六个辅助场之一，尚未重构完整 h_mu_nu。

6. 接入固定版本 pybhpt=1.0.0 的 s=±2 点粒子频域 Weyl 解。
   对公开表 rp=4/6/10/20 的两端振幅，时间积分后的最大相对误差约5.1e-15。
   特别注意：表的 Psi0/Psi4 数值与原始 pybhpt 曲率相差 i/omega。
   这被记录为显式的经验约定映射；原始曲率和时间积分分别保存，不能直接混用。
   文中定义的曲率与表的命名之间仍需在完整重构过程中逐步核对。
7. 背景求解扩展至论文的 alpha=.2；独立质量—荷检验及无径向节点检验通过。
   对应数值和背景图在 `alpha_02/`。

8. 自旋±1和紧支撑辅助标量的伴随点源积分已经实现。
   六个辅助场现在均可与公开幅度表比对；必须保留各自显式的归一化映射。
   实现、分部积分符号和独立验证见 [LORENZ_SOURCE_DERIVATION.md](LORENZ_SOURCE_DERIVATION.md)。
   这仍不等于完整 Lorenz 度规：张量组合的点源匹配及静态补全未完成。

9. 已实现非静态真空区域的三部分张量重构及 trace-driven kappa。
   `kappa_convergence.json`：对质量平方参数作四阶中心 resolvent 导数，
   改变步长、积分容差及远端边界后，场值与梯度的相对变化不超过4.3e-9。
   `metric_curvature_development.json`：独立坐标 Riemann 变分恢复输入
   psi0/psi4，相对差不超过6.4e-13；没有复用 GHP 重构算符作为检验。
   各部分通过真空 Einstein 和 Lorenz 检验，标量部分恢复预期 trace。
   进一步从自旋1规范向量独立计算 F=d xi，恢复输入的两个 Maxwell 标量，
   相对差小于3e-14。由此将错误的试验系数4修正为2。
   **完整度规仍缺点源匹配证明，不能接入生产环境源。**
10. 纠正中间复度规的验收定义：2406.12510v3 Sec. III A 明确先重构复解，
    再取实部。物理 Fourier 系数为 (h_m+conj(h_-m))/2。
    `metric_reality_development.json` 中原始复解的共轭差只是诊断，
    不应当作失败的现实性条件，也不能用取实部代替点源匹配验证。

当前环境及 Lorenz 模块共30项测试通过，包含新增的 Maxwell 反算检验。
背景图片是未受扰动云，不是论文尾迹图。

## Lorenz 重构公开资料调查

- [2306.16459v3](https://arxiv.org/abs/2306.16459)：作者注明 Mathematica 代码 available on request。
  用户确认没有该代码或度规数据。未联系作者。
- [2406.12510v3](https://arxiv.org/html/2406.12510v3)：公开源包含 `Operators.nb`、
  `amplitudes.dat`、正文 LaTeX；已下载到 ignored 的 `outputs/lorenz_reference/`。
  这些提供可读算符与独立数值基准，不是完整可直接运行的环境复现程序。
- [公开源包](https://arxiv.org/src/2406.12510v3) SHA256：
  `3bb7d110e55fb85a78aa56af85910ca0e1e5804c59016f504ec6ea913a0373a7`。
- 2406.12510v3 的非静态方案由六个带源 Teukolsky 场重构时间导数 Lie_T h。
  静态 m_g=0 仍须单独处理，不能通过除以零频率获得。
- Windows 有 WolframScript；当前 Python 实现未依赖运行外部 notebook。

下一阶段：核对完整复解的源约定，以点源连续性和导数跳跃
验证 Einstein 源，补齐静态和低多极度规。
之后接入已有协变源与 massive Green 求解器，推进上表图级验收。
目前缺作者程序不构成无法继续的理由：已有公开方程和数值表可支持独立实现。

## 重运行

```bash
.venv/bin/python src/environment_cloud.py
.venv/bin/python src/report_environment_background.py
.venv/bin/python src/report_environment_background.py --alpha .2
# report_lorenz_trace 需要上述公开源包位于 outputs/lorenz_reference/
.venv/bin/python src/report_lorenz_trace.py
.venv/bin/python src/report_lorenz_weyl.py
.venv/bin/python src/report_lorenz_spin1.py
.venv/bin/python src/report_lorenz_chi.py
.venv/bin/python src/report_corrector_adjoint.py
.venv/bin/python -m pytest tests/test_environment_cloud.py tests/test_environment_source.py tests/test_environment_radial.py tests/test_lorenz_trace.py tests/test_lorenz_weyl.py tests/test_environment_alpha02.py -q
```

所有当前数值结果只涉及此处明确列出的模块和参数，不代表全篇收敛或完整 EMRI 波形。
