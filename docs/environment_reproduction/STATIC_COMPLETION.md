# 静态Lorenz补全基：已构造，尚未匹配

依据2306.16459v3的completion表，实现A–G七个Kerr真空补全基。
它们还不是粒子诱导的物理静态解：系数、边界条件及与静态spin-2/trace
部分的联合匹配尚缺。约定M=1，纯规范项h_ab=2 nabla_(a xi_b)。

| 模式 | 构造 | 目标迹 |
|---|---|---|
| A | g_ab | 4 |
| B | 文献标量规范向量减2r_+² C | 6 |
| C | xi_a=(0,1/Delta,0,0) | 0 |
| D | xi^a=(t,a²cos²(theta)/Sigma,0,0)+nabla^a y | 2+2 Box y |
| E | partial_M g_ab - 2 nabla_a nabla_b y，固定a | -2 Box y |
| F | 文献时间线性旋转规范向量+nabla^a z | 4a |
| G | partial_a g_ab - 2 nabla_(a xi_b)，固定M | 0 |

D、F的向量显含t，但所得度规扰动静态。代码显式保留这些时间导数，
没有错误地用omega=0将它们置零。标量满足

    Box y = 2/(r_+-r_-) log((r-r_+)/(r-r_-)) = f(r),
    Box z = 2ar cos²(theta)/Sigma.

z使用文献闭式；y=y0(r)+y2(r)P2(cos theta)，满足

    (Delta y0')' = (r²+a²/3) f(r),
    (Delta y2')' - 6y2 = (2a²/3) f(r).

选择y0、y2及其一阶导数在参考半径6为零，得到明确特解。这未必与
文献按边界规则挑选的齐次自由度相同，不能直接套用文献的边界适配
系数。原文后续分离方程首行重复写z且遗漏log；这里从明确的Box y
方程乘Sigma推导径向式，并以独立张量检验验证。

在a=0、.6、.8771530275949366，r=3.5、9，theta=1.1，对全部七个
模式检查迹、Lorenz散度及独立线性化Einstein张量。最大绝对残差分别
约2.7e-15、1.4e-15、4.6e-15，见static_completion_validation.json。
这不证明物理边界条件或守恒荷归一化。

尚需：静态spin-2扇区、自由标量齐次项、轨道处联合匹配、质量/角动量
守恒荷独立核对和最终边界约定。本轮主项目测试58 passed。
