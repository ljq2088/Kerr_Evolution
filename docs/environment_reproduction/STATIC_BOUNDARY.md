# 静态 In/Up 边界与完整场候选

已实现 `environment_static_lorenz.StaticLorenzMode`，将静态 spin-2、
trace、匹配后的自由标量及补全项组合成协变 h_ab。当前是有限模态
候选；低阶源和边界的数值检查通过，不等于论文全部图的复现。

## 标量与二形式的 Green 解

M=1，b=sqrt(1-a²)，r+=1+b，x=(r-1)/b，lambda=j(j+1)。
标量径向算子及齐次基为

    L_j R = (Delta R')' - lambda R,
    P=P_j(x), Q=Q_j(x),
    W_s=Delta(P Q'-P'Q)=-b.

P 在视界正则，Q 在无穷远衰减。对真空分支上的扩展源 S(r)，

    R(r) = [Q(r) integral_(r+)^r P(s)S(s)ds
            +P(r) integral_r^infinity Q(s)S(s)ds]/W_s.

求一阶导数时积分端点项互相抵消，只对外面的 P/Q 求导。二阶导数
由 ODE 恢复。所有积分在粒子位置分段；这里的 S 是分支上的普通
函数，粒子处的分布项仍由完整度规匹配条件处理，不重复加入。

对二形式的 B_j：

    Delta B_j''-lambda B_j=S_B,
    U=Delta P_j', V=Delta Q_j',
    U'=lambda P_j, V'=lambda Q_j,
    W_B=U V'-U'V=lambda b.

将上一 Green 公式中的 P,Q,W,S 换成 U,V,W_B,S_B/Delta 即可。
j=0 二形式不存在。数值实现把两个基在 r0 归一化，避免高 j 时
原始 P/Q 幅值悬殊；Wronskian 作同样缩放。

代码 `lorenz_static_boundary.StaticGreen` 使用自适应复积分，真正
积分到无穷远，检验积分器状态。返回的积分误差估计不包括源误差、
Legendre 函数误差或近视界坐标变换中的浮点消减误差。
独立制造解 R=1/(r+.7) 和 B=Delta/(r+.7)^3（乘复常数）覆盖
Schwarzschild/Kerr、j=0/1/2/4、近视界与远区，恢复函数及两阶导数。

## 接入重构并保持局部匹配

`sourced_static_spin2(...,boundary='regular')` 的 chi 源使用两侧
各自的物理 Hertz 振幅；不能用同一个任意归一化 P/Q 分支跨轨道
构造源。`static_trace_metric(...,boundary='regular')` 同时为 kappa
和 B 指定上述边界。旧的 reference 选项仍保留，便于独立对照。

新旧扩展源相同，且两种 particular 的值及一阶导数都在 r0 连续，
故它们之差是跨轨道光滑的齐次规范场。实际检查 a=.6,r0=6,ell=2：
新旧差的度规值及一阶导数跳跃均小于7e-14。真空 Einstein 方程、
Lorenz 条件、trace 及极端 Weyl 曲率也经独立检验。

自由标量的两侧采用 P/Q。在 P(r0)=Q(r0)=1 的归一化下，令

    [kappa]=J0, [kappa']=J1,
    A_in=(J1-Q'(r0)J0)/(Q'(r0)-P'(r0)),
    A_up=A_in+J0.

由此构造 h_free,ab=nabla_a nabla_b kappa，既保留已求得的跳跃，
又固定边界齐次自由度。

补全中的 y2（乘 Legendre P2）也用 Green 解排除 r² 增长齐次项。
其源为 (2a²/3)f，f=b^-1 log((r-r+)/(r-r-))。参考基在 r0 取
y2=y2'=0，改基会改变 D/E，因此必须同步转换自由标量跳跃：

    delta J_(2,n) = -2(Delta c_D-Delta c_E)
                       y2^(n)(r0)/sqrt(5/(4pi)), n=0,1.

补全里的 y0 保留原参考约定。D+E 中 y 完全抵消，这一点另有测试。

## 从视界几何重推 F 补全条件

不能仅因为某组系数写在参考文中，就认为它适用于给定代码基底。
对这里按 2306.16459v3 表格与所印 z 实现的 F/G/C，直接采用文中的
C 系数多项式会在 advanced Kerr 的 h_rr 中留下二阶极点。
下述推导不通过拟合极点系数来确定答案。

记 rho_H²=r+²+a² cos²(theta)，Omega_H=a/(2r+)。
F 的时间依赖规范向量给出

    delta Omega_H[F]=2a Omega_H-1=-b,
    delta Omega_H[G]=partial_a Omega_H=1/(2b r+),
    delta Omega_H[C]=0.

因此 F+g G+c C 保持视界角速度要求 g=2b²r+。
在视界，所印 z 正则，Delta z' 的贡献为零；视界位置的变化为

    delta r_H[F]=2ab r+ cos²(theta)/rho_H²,
    delta r_H[G]=-a/b+a(r+ sin²(theta)+cos²(theta))/rho_H²,
    delta r_H[C]=-1/rho_H².

代入 g 后，角向依赖恰好消去，delta r_H=0 给出

    c=-2ab r+²,
    F_regular = F+2b²r+ G-2ab r+² C.

另一个正则组合是 D+E+aG+2r+²C。参考文所印 F 组合的 C 系数
a(r-³+r-²r++5r-r+²-7r+³)/6 与上述 c 相差 -4ab/3。
由于 C 的 advanced h_rr 主导项是 -1/[b(r-r+)²]，该差异留下
4a/[3(r-r+)²]。a=.6 时数值结果正是 .8/(r-r+)²。

改用几何推导的系数后，a=0,.6,.877 的测试均保持有界。这里明确
记录所印公式、当前基底与几何检验间的差异，不推断作者未公开代码
采用了何种修正或不同的齐次约定。

## 完整候选与当前边界证据

默认 Berndtson 选择令 Up 的 B/D/F 为零；In 系数取

    B_in=-Delta B, D_in=E_in=-Delta D, F_in=-Delta F,
    G_in=-a Delta D-2b²r+ Delta F,
    C_in=-2r+² Delta D+2ab r+² Delta F.

其他 Up 系数由已匹配的跳跃得到。它有非零内部质量/角动量补全，
不能当作零内部电荷的物理解。可选 physical 模式令 In 的 C/D/E/F/G
为零、B_in=-Delta B；这固定内部物理补全为零，但其 Lorenz 坐标
不渐近平坦。此选项尚未独立做完整边界与电荷积分验证。

以 a=.6,r0=6,L=8 的匹配数据构造默认候选：

* r-r+=.01,.001 时，advanced 坐标下最大协变分量为1.3494、1.3631，
  未出现之前的二阶极点；单独 BL 分量仍会随坐标奇性增大。
* r=30,100 时，以 L_a=(1,1,r,r) 缩放的最大 BL 分量为.05240、.01496，
  与远区衰减一致。两点不能替代解析渐近展开或全部边界误差研究。
* 全部五个四标架量的 j<=4 连续性残差2.58e-13，独立导数跳跃
  残差1.08e-11。原匹配系数没有用这些新边界场重新拟合。

完整数据在 `static_complete_boundary.json`，辅助场数据在
`static_auxiliary_boundary.json`。尚需论文截断范围、强场/近极端
参数、边界分辨率及全模态源的检验；不能由这些单个基准宣称整图复现。

主项目 tests 已通过70项，包括完整静态候选的真空检验、补全组合
的视界极点回归检验、自由标量跳跃及 Green 制造解检查。
