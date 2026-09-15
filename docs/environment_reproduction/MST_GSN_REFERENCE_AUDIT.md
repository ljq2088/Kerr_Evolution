# MST / GSN 适用性与本机独立参考审查

日期：2026-09-15。范围：只读检查参考文献、本机 GSN 代码与运行环境；本报告不声称已运行 GSN 对照，也不认定现有 32.66% 通量偏差已找到原因。

## 文献结论

“低频优先 MST、高频优先 GSN”是合理的数值策略，但并不是方程要求的、仅依赖一个固定 Mω 阈值的二分规则。Fujita–Tagoshi 已证明，允许复数重整化角动量 ν 后，MST 可超出早期实 ν 限制；实 ν 失效点依赖 s、ell、m、a。[Fujita–Tagoshi II](https://arxiv.org/abs/0904.3818)

Lo 的 GSN 工作将长程势变为短程势；GSN 在低频也可用，高频下 Riccati 相位积分往往更高效。该文同时展示：径向评估位置、渐近展开阶数、ell 和归一化都会影响比较；不能把 GSN 默认输出无条件当作真值。高 ell 默认边界展开的失败可以通过提高阶数修复。其 MST/GSN 对照应检查转换后的函数和 Wronskian，而非直接比较两套不同归一化的系数。[Lo 2024，尤其 III.3.2–III.3.3](https://arxiv.org/html/2306.16469v2)

ν 的求根/分支本身应独立检查；现代 monodromy 计算也有参数相关的相消限制，并非所有算法都天然稳定。[Nasipak 2025](https://arxiv.org/abs/2412.06503)

## 当前比较点的数值尺度（本地直接计算）

取 a=0.8771530275949366，r_p=20，M=1：

- Ω_p=0.011071760582072043。
- r_+=1.4802109600800843，Ω_H=0.29629324847975713。
- 视界表面引力 κ_H=0.1622103109053129。
- 度规驱动频率 ω_g=m_g Ω_p；m_g=1,…,6 对应 0.01107176058,…,0.06643056349。
- 视界频率 p_g=ω_g−m_g Ω_H=−0.2852214878976851 m_g。
- Frobenius 视界相位指数的系数 p_g/(2κ_H)=−0.8791718797215596 m_g。
- MST ε=2Mω_g=0.022143521164144087 m_g，而 τ_MST=(ε−m_g a)/sqrt(1−a²)=−1.7804872806072636 m_g。

因此度规 ω_g 属于低频，但不能据此认为 Kerr 视界附近也是小频率/小相位问题。m_g 增大时 τ_MST 与视界相位增长。应保留完整的旋转边界指数，检查导数、正负自旋变换、端点展开误差以及实际源积分节点处的函数精度。

环境中的有质量标量频率是 ω_c+m_g Ω_p；其远区波数为 sqrt(ω²−μ²)，闭合通道则为 sqrt(μ²−ω²)。这与无质量的度规 Teukolsky 函数是不同径向问题。本机 GSN 参考接口未提供 μ 参数，不能用它直接替代 massive scalar Green 函数，也不能把云的约0.296频率误当作度规驱动频率。

## 本机参考代码

位置：`/home/ljq/code/GeneralizedSasakiNakamura.jl`。
版本：`0.7.1`。
提交：`a403cf51e3d1ca85f4923f011b8a131b14d14fd5`。
官方仓库：<https://github.com/ricokaloklo/GeneralizedSasakiNakamura.jl>。

这是比当前官网新版本早的本机代码。读取本机代码所得：

- 支持 s=0,±1,±2；质量 M=1。
- 实频 `method="auto"` 选 `Riccati`；`linear` 为独立的同一 GSN 方程积分形式。
- 默认 ODE 为 Vern9，ComplexF64，绝对/相对 tolerance=1e-12。
- 默认数值积分 r_*∈[−50,1000]；地平线级数3阶、无穷远6阶。
- 高层成对接口会增加边界展开阶数，并以渐近入射系数比的恒等式判断；该判断在 `Solutions.jl:744` 的默认容差为绝对1e-6，不能直接当作1e-12的全局误差保证。
- 本机 `GeneralizedSasakiNakamura.jl:582–584` 在所需边界阶数超过50时打印警告后 `continue`，而未改变条件或退出：若一直不满足 sanity，可能形成无限循环。未修改参考库。独立审查可避开此高层循环，调用显式边界接口并自行检查。
- `Teukolsky_radial` 会应用 GSN→Teukolsky 振幅转换，并返回单位 Teukolsky transmission 归一化。R(r)与dR/dr均保存在 `.Teukolsky_solution(r)` 中；只调用 `Rin(r)` 返回第一个分量。

可复用的显式接口（尚未执行）：

```julia
using GeneralizedSasakiNakamura
s = -2; l = 2; m = 1
a = 0.8771530275949366
omega = m / (20.0^1.5 + a)
Rin = Teukolsky_radial(s, l, m, a, omega, IN, -50.0, 1000.0;
    method="linear", horizon_expansion_order=15,
    infinity_expansion_order=20, tolerance=1e-12)
Rup = Teukolsky_radial(s, l, m, a, omega, UP, -50.0, 1000.0;
    method="linear", horizon_expansion_order=15,
    infinity_expansion_order=20, tolerance=1e-12)
Rin_value, Rin_prime = Rin.Teukolsky_solution(20.0)
Rup_value, Rup_prime = Rup.Teukolsky_solution(20.0)
```

上述阶数与范围是起始诊断设置，并非已证明收敛的选择；还需逐一改变展开阶数、端点、ODE 容差，并与 Riccati 结果比较。

本轮检查 Windows/WSL PATH、`~/.julia`、`~/.juliaup`、`/opt`、`/usr/local/bin` 没有发现 Julia 可执行或现成 Julia 环境。conda 的包缓存中有 `pyjuliacall`、`pyjuliapkg` 记录，但这不是 Julia 运行时已安装的证据；未通过导入自动下载任何依赖。没有安装依赖，也没有修改本机参考仓库。

## 推荐的误差审查判据（本报告提出的实现审查方案）

1. 对实际生产 s、ell、m、a、ω 和整段源积分半径，分别强制 MST 与另一独立径向方法；比较归一化无关的 R'/R 和固定半径归一后的复值曲线。
2. 检查 Δ^(s+1)(Rin Rup'−Rup Rin') 沿半径的相对漂移；归一化一致后再比较其绝对值与两端散射恒等式。
3. 改变边界位置和渐近展开阶数；只有固定 ODE tolerance 的单次测试不能排除初值截断误差。
4. 检查解析 radial ODE 递推导数与正负自旋的 Teukolsky–Starobinsky 恒等式。ODE内部根据同一方程产生的二阶导数残差并不构成独立精度验证。
5. massive scalar 单独检查质量阈值、衰减/振荡分支、远区级数与 Wronskian 相消；不能以 massless GSN 验证取代。
6. 最后将径向差异传播到复源、ZH与通量。若只看局部函数相对误差，不足以证明或排除32.66%的最终通量差异。
