# Lorenz 辅助场源与伴随计算

参考固定到 [2406.12510v3](https://arxiv.org/html/2406.12510v3)，原始源包哈希见已有验证记录。
本模块提供六个辅助场的渐近幅度输入；完整度规、静态补全和文章尾迹仍未完成。

## 实现和独立校验

`lorenz_jet.py` 用二元 Taylor 系数表示 r、theta 导数，避免对高阶算符逐层有限差分。
`lorenz_ghp.py` 从 Kinnersley 四标架的协变导数直接计算自旋系数，显式跟踪 GHP 权重。
`lorenz_mode_jet.py` 从径向/角向 Teukolsky ODE 递推任意所需的局部导数。

独立检查包括：四标架内积、自旋系数、O_0=-Box/2、五种自旋分离方程及其一阶导数，
以及 conformal Killing-Yano 恒等式。使用的 CKY 恒等式为

    nabla_c f_ab = g_bc T_a - g_ac T_b.

这保持 a,b 反对称。不能照搬源 TeX 中第二项 g_ab T_c 的指标排列。
spin=-1 的 GHP 角算符使用与 priming 对称及分离方程一致的符号；独立 ODE 检查通过。

`lorenz_corrector.py` 逐项实现附录 B 的 N 修正算符，**不吸收 8pi**。
由 trace-free 和 N^dag=-N 可知 N[g f]=0；该测试通过。
另用相反 Fourier 频率的紧支撑对称张量做双线性积分（不取复共轭），验证

    integral [V^{ab} N(U)_ab + U^{ab} N(V)_ab] sqrt(-g) = 0.

Gauss 阶数8至12，归一化缺陷从1.09e-6降至1.06e-11，见 `corrector_adjoint.json`。

## 点粒子源的伴随积分

令质量比已剥除的四速度为 u，粒子 Fourier 源 T^{ab}=u^a u^b/(u^t Sigma sin theta) delta。
对任意伴随测试协向量 V_a，定义 D_ab=nabla_b V_a，U_ab=D_(ab)。
从 j_SE=8pi div N[T]-j_trace 分部积分得到

    integral V^a j_SE_a sqrt(-g)
      = (8pi/u^t) [u^a u^b N(U)_ab - f^{ab} D_ab - V_t] at particle.

这里因 N^dag=-N，修正项取正；trace T 沿粒子的收缩为 -1。
`current_adjoint_contraction` 实现方括号，调用者补 8pi/u^t 和 Wronskian。

对 spin=s 的径向算符，Wronskian 为 Delta^(s+1)(R_In R_Up'-R_Up R_In')。
伴随 Green 核（不含 Wronskian）为

    -2 Delta^s zeta^(|s|-s) R_other S_s(theta) exp(+i omega t-im phi).

-2 来自 GHP Teukolsky 算符与分离径向方程的主部关系。
角函数与径向数据仍取原始 (s,m,omega)，只有作为测试函数的 Fourier 相位反向。
`lorenz_spin1.py` 将该核送入 Maxwell S^dag，再使用上述点粒子积分。
源中的局域接触项不贡献两端幅度；完整度规在粒子处的分布抵消仍需另行实现。

## 紧支撑辅助标量

计算的是 chi_compact=chi_DKWSE+chi_SE；完整度规还需要 kappa，满足 Box kappa=Lie_T h。
由于 h 在轨道两侧延伸，kappa 不是普通紧支撑点源的齐次 In/Up 拼接问题。

从展开方程出发：

    Box chi_DKWSE = div xi_DKW,
    Box chi_SE = div(f J_SE),
    J_SE = Lie_T^-1 (j_SE-j_DKWSE).

分部积分 div(f J) 给出 **-i/omega** 乘以 j 的伴随收缩。
必须将源中的 chi^dag 移到 Green 核上后改成 chi，而不是再次使用 chi^dag。
`lorenz_chi.py` 直接构造展开的 j_DKWSE 伴随作用，保留全部三项；
没有用一个拟合常数修正原先未通过的压缩公式实现。

## 与公开表的显式映射

公开 `amplitudes.dat` 中的命名不能直接当作本实现的物理场定义。四个轨道
rp=4,6,10,20、a=.6、ell=m=2 的比对得到以下统一关系：

| 表中列 | 本实现的原始幅度到表的映射 | 最大相对差 |
|---|---|---|
| h | 原值 | 约2.4e-9 |
| Psi0, Psi4 | i/omega | 约5.1e-15 |
| Phi0, Phi2 | sqrt(2)/omega² | 约1.2e-13 |
| Chi | -i/(2 omega) | 见 `lorenz_chi_validation.json` |

映射只用于公开表的比对；原始幅度单独保存。完整重构应从明确的微分方程和
Fourier 约定组装，并用 Lorenz 条件、Einstein 方程和现实性关系继续核验，
不能把表的重标度再重复施加到物理曲率或规范向量。

## 下一步

1. 从原始 Weyl 场构造 AAB 与 DKW 规范向量；在真空区域核对 spin-2 度规部分。
2. 用原始 spin-1 和 chi_compact 构造其余规范部分，求 trace-driven kappa。
3. 验证完整 Fourier 度规的 Lorenz 条件、真空 Einstein 残差、轨道连续性/源跳跃、
   h_m=conjugate(h_-m)；静态 m=0 和低多极补全单独处理。
4. 此后才接入已验证的云 Hessian 和 massive Green 模块，执行全文章图级复现。
