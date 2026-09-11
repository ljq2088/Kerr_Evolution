# 完整复现目标与当前证据

目标保持为完整复现 2501.09806v1 正文及补充材料的数值结果。
目前状态：进行中，不能由阈值一致或下列模块测试推断论文已复现。

最新图间核对发现 Fig.6 的 ell<=6 模态和除以 alpha^6，与 Fig.2
同半径总通量相差约4.2e-5，支持两图数值采用不同归一化的推断。
详见 `PAPER_NORMALIZATION.md`。scalar(3,3) 剩余约5%的差异尚未
解决；云质量测度的独立排查不足以解释它。所有图级目标仍未验收。

近阈值实际源试算发现固定远区边界的逆r级数失效，已加入Whittaker–Coulomb
边界并用独立远边界计算交叉核对。见 `THRESHOLD_BOUNDARY.md`；旧失效
通量明确保留为历史诊断，不能用于图5验收。

Schwarzschild复频率背景、KS切片能量归一化及冻结时间衰减的KG缺陷
已独立实现与检查，见 `SCHWARZSCHILD_BACKGROUND.md`。Schwarzschild
后续已接入显式标记的Schwarzschild冻结背景受迫通道，完成首轮有限
分辨率积分与球对称选择规则检验；尚未完成收敛与总通量，图2/3未完成。
静态A–G补全基已通过真空身份检验，详见 `STATIC_COMPLETION.md`，
但匹配系数与物理边界条件尚未确定，不能视为静态扇区已经完成。

静态spin-2 Hertz–Lorenz分支及点粒子曲率归一化已实现并检查Psi0/Psi4，
见 `STATIC_SPIN2.md`。联合trace、自由标量和补全基的轨道匹配仍未完成。
已实现联合局部匹配并修正圆轨道等距对称性。采用论文的四标架
条件后，Kerr低阶独立导数检验呈二阶外推收敛，误差约5e-11；
图1强场参数也通过低阶局部检验。模态截断边缘及物理边界尚未完成，
详见 `STATIC_MATCHING.md`，不能据此宣布任一图已复现。
现已接入辅助场的 In/Up Green 解及自由标量边界，并组合成静态场
候选。默认 Berndtson 基准通过低阶源和有限边界检查；F 补全条件
经视界几何重推以消除实际极点，详见 `STATIC_BOUNDARY.md`。
标量(0,0)主导视界通道已完成真实源积分与首轮独立加密，
近视界截断仍限制精度，详见 `SCALAR00_RESPONSE.md`。
已有视界渐近积分补偿及阶数/区间交叉检验，见 `HORIZON_TAIL.md`；
它减少了单通道的内截断误差，但其余收敛和整图目标仍未完成。

## 图级验收

| 目标 | 所需证据 | 当前状态 |
|---|---|---|
| 图1，alpha=.3，rp=3.5，赤道与方位切面尾迹，标量l>=2 | 完整 Lorenz 点粒子源、受迫场、坐标约定、收敛图 | 未完成 |
| 图2，alpha=.2/.3，Kerr 与 Schwarzschild 两端通量 | 正确云归一化，无穷远标量l<=6、视界标量l<=5，逐轨道输出与边界收敛 | 未完成 |
| 图3，相对通量误差及旧 Schwarzschild 实现比较 | 同一归一化及规范处理，原始曲线/独立求解器核对 | 未完成 |
| 图4，标量/GW 通量比 | q=1e-6、eta=.1，与作者确认或解释 v1 epsilon 幂次差异 | 未完成 |
| 图5，rp=41.6/41.8 两幅尾迹 | 独立求解两轨道，radiative/bound 过渡，无伪造源 | 阈值位置已核对；场图未完成 |
| 图6，rp=20 的逐模通量 | 模态奇偶与高 l 收敛、主/次模幅度 | 未完成 |
| 图7，rp=20 粒子处 l^-2 收敛 | 点粒子真实源、场连续性、高 l 模和 | 未完成 |
| 共振解释与位置 | 束缚谱和有源响应，普通/轨道有效能流区分 | 未完成 |

## 已实现且验证的模块

论文补充材料 S.5 的度规求和截断是 ellmax=18；图1/5的 l>=2 是
绘制的标量场截断，两者不能混淆。当前度规 ell<=4 的单通道试算
不满足论文的完整截断要求，仍需加密和模态求和。

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

11. 显式加入同一复电流的反自对偶 Maxwell 重构，保留不能在该通道中
    丢弃的 DKW 源积分。初步 h_tt 角向投影的连续性和导数跳跃显著改善，
    推导和数值边界见 `LORENZ_METRIC_DEVELOPMENT.md`。这尚不构成所有
    张量分量、各参数或静态度规的点源验收。

12. 全部10个独立张量分量的两类角向投影正在验证；Taylor 阶数必须
    至少为8才能可信地比较全部分量的二阶导数。ellmax=6 仍有明显截断误差，
    不能仅根据 h_tt 的良好结果宣称完整点源匹配通过。
    零频 m=0 trace 已以 Legendre Green 函数独立实现，静态完整度规尚未完成。

当前环境及 Lorenz 模块共40项测试通过，包括静态 trace 规范部分及非静态偶极检验。

13. 已连接真实重构度规与云的协变源，不再仅用测试张量验证接口。
    `actual_lorenz_source_development.json`：alpha=.3、rp=20、r=10，
    m_g=2、度规 ell<=4、标量(ell,m)=(3,3)，6点与10点角积分相差约2.3e-9。
    这是一个有限模态的局部有源方程输入，不是受迫场、通量或论文图的复现。
14. 一个真实非静态通道已完成两套径向源积分，并输出有限源区间的
    受迫响应、两端幅度和通量。每面板2点至4点仍有明显通量变化，尚未收敛。
    详情及归一化见 `FORCED_MODE_DEVELOPMENT.md`；图级验收状态仍未完成。
15. 静态 trace 对应的真空 Lorenz 规范部分已按2306.16459独立构造，
    仍缺自旋2、自由齐次规范项及质量/角动量补全的匹配。
    非静态 ell=1 的标量/矢量部分已接入，用于 m_g=1 通道；该通道
    的完整点源匹配和标量通量尚未验收。
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
