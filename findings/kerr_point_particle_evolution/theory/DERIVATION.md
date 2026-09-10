# W02.1 点粒子有源 Teukolsky 方程：阶段 1.1--1.6

- 任务：T006
- 状态：owner 已接受当前连续理论/符号检查范围；A1 与 boundary-adjacent finite-difference closure 仍开放，故不作完整 Gate A/production claim
- 范围：从未投影 $(T,R,\theta,\Phi)$ 方程到 verified interior smooth source，并给出 A6 final PDE/first-order/axis/SCRI theory contract draft

## 推导路线

本文遵循唯一的推导链：

```text
1. 固定 Ripley 真空算符和演化场
2. 写出 Sasaki BL/Kinnersley 有源方程
3. 用四旋标权重连接两种场变量
4. 把 BL 微分算符拉回到双曲型紧致化坐标
5. 用 Sasaki 左端复现 Ripley 左端，确定整体比例
6. 对 Sasaki 右端施加同一拉回，得到四块源
7. 在常 T 超曲面上构造点粒子应力-能量张量
8. 得到三个 Kinnersley 投影和运动 delta 函数 Jacobian
9. 合并为未投影的候选有源方程
10. Fourier 投影到固定 m，并变换 y=-cos(theta)
11. 冻结 fixed-m、(T,R,y) exact distributional equation
12. 以 dR dy measure 引入 Gaussian，并指定 radial off-worldline extension
13. 解析计算 worldline、amplitude、phase 和 moving-Gaussian time jets
14. 在四个 source blocks 中消去 coordinate partial_T
15. 固定 post-ISCO initial-data 与唯一 hyperboloidal trajectory interface
16. 把 sourced second-order PDE 化为连续 P/Q first-order system
17. 推导 characteristics、axis regularity 和 SCRI+ projection contract
```

正文只保留这条推导所需的关键等式。完整 Ripley 算符、Sasaki 四块源和点粒子分量分别见 [附录 A](appendix/RIPLEY_OPERATOR.md)、[附录 B](appendix/SASAKI_SOURCE.md)和 [附录 C](appendix/POINT_PARTICLE_PROJECTIONS.md)；fixed-$m$/$y$ 与 Gaussian/time-jets 的长推导见 [附录 D](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md)和 [附录 E](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md)。约定汇总见 [CONVENTIONS.md](CONVENTIONS.md)。

## I. Ripley 目标方程

采用双曲型紧致化坐标

$$
T=t+r_*-2r-4M\ln(r/M),
\qquad
R=\frac{L^2}{r},
\qquad
\Phi=\phi+q(r),
\quad q'(r)=\frac{a}{\Delta}.
$$

$R=0$ 对应未来零无穷远。由坐标定义的渐近关系直接有

$$
T\big|_{\mathscr I^+}=u=t-r_*.
$$

Ripley 的正则四旋标场为

$$
\Psi_4^{\rm R}=R\psi_4,
$$

其未投影真空方程记为

$$
\mathcal O_{\rm R}[\psi_4]=0.
$$

$\mathcal O_{\rm R}$ 的完整展开见 [附录 A](appendix/RIPLEY_OPERATOR.md)。当前 Mathematica 比较把该展开作为目标输入；它尚未独立完成从 Ripley Eq. (22) PDF 到表达式的转录核验。

## II. 从 Sasaki 方程确定有源项

### A. 原始有源方程

Sasaki--Tagoshi 的 Kinnersley 四旋标主场满足

$$
\phi_{\rm ST}=\rho^{-4}\Psi_4^{\rm K},
\qquad
\mathcal O_{\rm ST}[\phi_{\rm ST}]
=4\pi\Sigma\widehat T_{\rm ST},
$$

其中

$$
\widehat T_{\rm ST}=2(B'_2+B_2^{*\prime}),
$$

其源的层级关系为

$$
T_{\mu\nu}
\longrightarrow
\left\{T_{nn},T_{\bar mn},T_{\bar m\bar m}\right\}
\longrightarrow
\left\{B'_2,B_2^{*\prime}\right\}
\longrightarrow
\widehat T_{\rm ST}.
$$

$B'_2$ 和 $B_2^{*\prime}$ 是三个应力-能量张量投影的有序微分组合，并非新的物质变量。四个源块的系数和嵌套次序见 [附录 B](appendix/SASAKI_SOURCE.md#四个有序-source-blocks)。

### B. 场变量与坐标拉回

由两套四旋标之间的 boost/spin 权重，

$$
\phi_{\rm ST}
=\frac{\Delta^2}{4}\Psi_4^{\rm R}
=\frac{\Delta^2R}{4}\psi_4.
$$

令

$$
h'(r)=\frac{r^2+a^2}{\Delta}-2-\frac{4M}{r}.
$$

则 BL 微分算符的拉回为

$$
\partial_t=\partial_T,
\qquad
\partial_\phi=\partial_\Phi,
$$

$$
\partial_r\big|_{t,\phi}
=h'(r)\partial_T
-\frac{L^2}{r^2}\partial_R
+\frac{a}{\Delta}\partial_\Phi.
$$

把该拉回作用于 Sasaki 左端，得到与 Ripley 算符的恒等关系

$$
\boxed{
\mathcal O_{\rm ST}
\left[\frac{\Delta^2R}{4}\psi_4\right]
=-\frac{\Delta^2R}{4}\mathcal O_{\rm R}[\psi_4]
}.
$$

因此，无需在中间表达式层面强行统一 Ripley 与 Sasaki 的 conventions；只要从 Sasaki 有源方程出发一致地变换左右两端，并证明左端等于实际采用的 Ripley 数值算符，右端就给出与该算符配套的 source，但当前 A1 的 signature/unit bridge 仍需闭合。

**MCP 核验 A1/A2（场与左端算符）。** Mathematica 中分别输入 Sasaki 算符和 Ripley 目标算符，再对 Sasaki 场代入 $\Delta^2R\psi_4/4$ 并施加上述坐标拉回。比较的整式 residual 为

$$
\mathcal R_{\rm op}[\psi_4]
=\mathcal O_{\rm ST}
\left[\frac{\Delta^2R}{4}\psi_4\right]
+\frac{\Delta^2R}{4}\mathcal O_{\rm R}[\psi_4].
$$

由于这是线性二阶算符，程序用算符中实际出现的 12 个独立 jets 隔离各项系数：

$$
\mathbb J=
\left\{
\psi_4,
\psi_{4,T},\psi_{4,R},\psi_{4,\theta},\psi_{4,\Phi},
\psi_{4,TT},\psi_{4,TR},\psi_{4,RR},
\psi_{4,\theta\theta},\psi_{4,\Phi\Phi},
\psi_{4,T\Phi},\psi_{4,R\Phi}
\right\}.
$$

这些 jets 由 12 个 polynomial test fields 在一般点逐一激活。`OperatorResidualsByJet` 的 12 项全部严格化简为 $0$；坐标 Jacobian、二次径向 chain rule 和 $\widehat{\mathcal J}_+$ 的系数 residual 也全部为 $0$。场因子 `A1MasterFieldResidual` 同样为 $0$，但它以已输入的 boost/spin relation 为前提。执行证据见 [A1/A2 结果](../../../data/kerr_point_particle_verification/results/A1_A2_operator.json)。

### C. 对右端施加同一拉回

Sasaki 源中的两个基本微分算符变为

$$
\widehat{\mathcal L}_s
=\partial_\theta-i\csc\theta\,\partial_\Phi
-ia\sin\theta\,\partial_T+s\cot\theta,
$$

$$
\widehat{\mathcal J}_+
=-\left(2+\frac{4M}{r}\right)\partial_T
-\frac{L^2}{r^2}\partial_R.
$$

$\widehat{\mathcal J}_+$ 中的 $\partial_\Phi$ 在固定 $m$ 之前已经相消。把 $\mathcal L_s$、$\mathcal J_+$ 分别替换为 $\widehat{\mathcal L}_s$、$\widehat{\mathcal J}_+$，并保持四个源块的原有次序，得到作用于三个 Kinnersley projections 的线性 source functional：

$$
\boxed{
\begin{aligned}
\widehat T_H
\left[T_{nn}^H,T_{\bar mn}^H,T_{\bar m\bar m}^H\right]
={}&
\widehat{\mathcal B}_{nn}^{(LL)}[T_{nn}^H]
+\widehat{\mathcal B}_{\bar mn}^{(LJ)}[T_{\bar mn}^H]\\
&+\widehat{\mathcal B}_{\bar m\bar m}^{(JJ)}[T_{\bar m\bar m}^H]
+\widehat{\mathcal B}_{\bar mn}^{(JL)}[T_{\bar mn}^H].
\end{aligned}
}
$$

上标 $(LL),(LJ),(JJ),(JL)$ 按外层到内层标记四块的 operator composition；每个 $\widehat{\mathcal B}$ 的完整 prefactor 和 nesting 见 [附录 B](appendix/SASAKI_SOURCE.md#拉回后的-source)。因此 $\widehat T_H$ 不是应力-能量张量的一个坐标分量，而是由其三个 tetrad projections 构造出的 Teukolsky source。对任意 matter source，

$$
-\frac{\Delta^2R}{4}\mathcal O_{\rm R}[\psi_4]
=4\pi\Sigma\widehat T_H,
$$

即

$$
\boxed{
\mathcal O_{\rm R}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}\widehat T_H
}.
$$

**MCP 核验 A2（右端 source）。** Mathematica 从原始 Sasaki BL four-block source 出发，把三个 projections 作为任意 composed functions，直接执行 coordinate composition、chain rule 和 product rule；生成路线不使用 hatted operators、理论 derivative map 或 target block template。所得完整 hyperboloidal source 与单独录入的 theory $\widehat T_H$ 在四块、两组 $B$ 和总 endpoint 上 residual 全部为 $0$，保存文件与 root replay 均成功。该结果闭合了 theory BL start 到 theory hyperboloidal endpoint 的代数；primary PDF transcription 以及 A1 signature/sourced-unit bridge 仍开放。执行证据见 [A2 结果](../../../data/kerr_point_particle_verification/results/A2_four_source_blocks.json)和 [左端 operator 结果](../../../data/kerr_point_particle_verification/results/A1_A2_operator.json)。

## III. 赤道 Kerr 点粒子源

本阶段继承 Section II 已给出的 A2 source functional；不重新推导其外部四个 operators。给定 W02.3 将提供的 $(E,L_z,r_0)$，定义

$$
B=L_z-aE,
\qquad
P(r)=E(r^2+a^2)-aL_z,
$$

$$
C(r)=r^2+B^2,
\qquad
\mathcal R(r)=P(r)^2-\Delta(r)C(r).
$$

`P08` Eqs. (27)--(28) 与 `P09` Eqs. (2.1)--(2.3) 给出 Kerr geodesic first integrals。对 $Q=0$ 的赤道 future-directed inward branch，

$$
\boxed{
u^r=-\frac{\sqrt{\mathcal R}}{r^2},
\qquad u^\theta=0
},
$$

$$
\boxed{
u^t=\frac1{r^2}
\left[aB+\frac{(r^2+a^2)P}{\Delta}\right],
\qquad
u^\phi=\frac1{r^2}
\left[B+\frac{aP}{\Delta}\right]
}.
$$

这些 expressions 满足 $g_{\mu\nu}u^\mu u^\nu=-1$、$E=-u_t$ 和 $L_z=u_\phi$。取 $\mathcal R\ge0$、$u^t>0$；精确 turning point 与 horizon 分别作为 branch/domain boundary 记录。

坐标 Jacobian 给出

$$
u^T=u^t+h'u^r,
\qquad
u^R=-\frac{L^2}{r^2}u^r,
\qquad
u^\Phi=u^\phi+\frac a\Delta u^r,
\qquad u^y=0.
$$

由

$$
\boxed{
K=\frac{P-\sqrt{\mathcal R}}{\Delta}
=\frac{C}{P+\sqrt{\mathcal R}}
}
$$

得到适合 future horizon 的 regular form

$$
\boxed{
u^T=\frac1{r^2}
\left[
aB+(r^2+a^2)K
+\left(2+\frac{4M}{r}\right)\sqrt{\mathcal R}
\right]
},
$$

$$
\boxed{
u^R=\frac{L^2\sqrt{\mathcal R}}{r^4},
\qquad
u^\Phi=\frac{B+aK}{r^2},
\qquad u^y=0
}.
$$

这里要求 $P+\sqrt{\mathcal R}\ne0$；future-horizon crossing 取 $P(r_+)>0$，此时 $K\to C(r_+)/[2P(r_+)]$。

给定 $r_p(T_0)=r_0$ 后，worldline 由

$$
\frac{dr_p}{dT}=\frac{u_p^r}{u_p^T},
\qquad
\frac{dR_p}{dT}=\frac{u_p^R}{u_p^T},
\qquad
\frac{d\Phi_p}{dT}=\frac{u_p^\Phi}{u_p^T}
$$

确定；$E,L_z,r_0$ 的具体 ISCO/transition prescription 留给 W02.3。

`P08` Eq. (26) 给出这一 tensor 的 BL time-sliced form，Eqs. (29)--(32) 给出其 Kinnersley projections。其 invariant starting point 是

$$
T^{AB}(X)=\mu\int\frac{u^Au^B}{\sqrt{-g_H}}
\delta^{(4)}\!\left(X-z(\tau)\right)d\tau.
$$

在 worldline 与 constant-$T$ slice 唯一且 transverse 相交时，$\delta(T-T_p(\tau))$ integration 给出 $1/|u_p^T|$。采用 $u_p^T>0$ 的 future branch，并用

$$
\sqrt{-g_H}=\frac{\Sigma r^2\sin\theta}{L^2},
$$

得到

$$
\boxed{
T_H^{AB}
=\frac{\mu u_p^Au_p^B}{\sqrt{-g_H}\,u_p^T}
\delta(R-R_p)\delta(\theta-\theta_p)\delta(\Phi-\Phi_p)
}.
$$

三个 deltas 相对于 $dR\,d\theta\,d\Phi$ 定义。从 BL time-sliced tensor 独立拉回时，fixed-$T$ root derivative 为 $u_p^T/u_p^t$，并有

$$
\frac1{u_p^t}\delta\!\left(r-r_p(t)\right)
=\frac{L^2}{r_p^2u_p^T}
\delta\!\left(R-R_p(T)\right),
$$

而 azimuthal delta 在 radial support 上变为 $\delta(\Phi-\Phi_p(T))$，不产生额外 Jacobian。

定义 Kinnersley contractions

$$
N_p=n_A^{\rm K}u_p^A,
\qquad
\bar M_p=\bar m_A^{\rm K}u_p^A,
$$

以及

$$
\mathcal W_p^{(\theta)}
=\frac{\mu L^2}
{\Sigma_pr_p^2\sin\theta_p\,u_p^T}
\delta(R-R_p)\delta(\theta-\theta_p)\delta(\Phi-\Phi_p).
$$

则

$$
\boxed{
T_{nn}^{H,\rm pp}=\mathcal W_p^{(\theta)}N_p^2,
\qquad
T_{\bar mn}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_pN_p,
\qquad
T_{\bar m\bar m}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_p^2
}.
$$

在赤道 inward branch 上，

$$
\boxed{
N_p=-\frac{P_p-\sqrt{\mathcal R_p}}{2r_p^2}
=-\frac{\Delta_pK_p}{2r_p^2},
\qquad
\bar M_p=\frac{i(aE-L_z)}{\sqrt2\,r_p}
}.
$$

把 $T_{nn}^{H,\rm pp},T_{\bar mn}^{H,\rm pp},T_{\bar m\bar m}^{H,\rm pp}$ 逐项代入 A2 四块，定义具体的 $\widehat T_{H,A}^{\rm pp}$，得到

$$
\boxed{
\widehat T_H^{\rm pp}=\sum_{A=1}^4\widehat T_{H,A}^{\rm pp},
\qquad
S_{\rm R}^{\rm pp}
=-\frac{16\pi\Sigma}{\Delta^2R}
\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}
},
$$

$$
\boxed{
\mathcal O_{\rm R}[\psi_4]=S_{\rm R}^{\rm pp}
}.
$$

完整 Kinnersley covectors、contractions、absolute-value/domain assumptions 以及四个已经代入 $\mathcal W_p^{(\theta)}N_p^2$、$\mathcal W_p^{(\theta)}\bar M_pN_p$、$\mathcal W_p^{(\theta)}\bar M_p^2$ 的 nested blocks 见 [附录 C](appendix/POINT_PARTICLE_PROJECTIONS.md)。

**MCP 核验 A3。** Mathematica 分为 geodesic velocities、invariant constant-$T$ route、BL delta pullback route 与 tetrad projections 四个独立 stages。BL/hyperboloidal four-velocity、timelike normalization、$K$ identity、future-horizon regular limit、determinant、absolute-value branches、全部 spatial deltas、covectors、$N_p,\bar M_p$ 及三条完整 projection route residuals 均为 $0$；summary 给出 `A3FullEndpointVerified=True`。A3 在 $T_{nn}^{H,\rm pp},T_{\bar mn}^{H,\rm pp},T_{\bar m\bar m}^{H,\rm pp}$ 处闭合，外部 four-block functional 由 A2 覆盖而不重复验证。结果见 [A3 endpoint summary](../../../data/kerr_point_particle_verification/results/a3/a3_endpoint_summary.json)。该结论是 local algebraic closure：global worldline monotonicity 与 global unique slice intersection 未被验证，并继续继承 A1 convention/tetrad boundary。

## IV. 代数核验与开放边界

MCP 是 Mathematica Kernel 的执行通道，不是另一套代数方法。现有 A1--A6 artifacts 的共同证据结构为

```text
物理起点、conventions 与 assumptions
  -> Kernel 独立生成计算终点
  -> 单独录入理论终点
  -> 选择适合当前对象的精确表示与化简
  -> 构造 residual = computed endpoint - theory endpoint
  -> 要求所声明范围内的 residual 精确等于 0
  -> MCP 返回结果并保存为 JSON
```

对任意待验证关系 $F=G$，实际判据可概括为

$$
\mathcal R
\equiv
\operatorname{ExactReduce}_{\mathcal A}[F-G],
\qquad
\text{pass}\iff\mathcal R\equiv0,
$$

其中 $\mathcal A$ 是随结果明确记录的实数域、非零分母、branch 和坐标范围假设；$\operatorname{ExactReduce}$ 不是固定程序函数，而表示按对象选择可审阅的精确代数路线。A1--A3 主要使用 Together/FullSimplify；A4 的大分布源使用完整求导后的 Expand/jet coefficients；A5 使用 common-jets 与四块 formal-spatial-derivative interface 的模块化比较；A6 只做一阶化、characteristics 和 axis/projection 的 targeted checks。这里的 $0$ 是精确零，不是抽样点上的数值近零。核验没有使用 PowerExpand，transport/debug failures 与数学 residual 分开记录。

代码与原始结果集中保存在 `code/kerr_point_particle_verification/checks/` 与 `data/kerr_point_particle_verification/results/`：

| 核验对象 | 自包含输入 | MCP 返回数据 |
|---|---|---|
| A1/A2：场因子、坐标变换、完整左端算符、整体 source factor | [A1_A2_operator.wl](../../../code/kerr_point_particle_verification/checks/A1_A2_operator.wl) | [A1_A2_operator.json](../../../data/kerr_point_particle_verification/results/A1_A2_operator.json) |
| A2：四个有序 source blocks | [A2_four_source_blocks.wl](../../../code/kerr_point_particle_verification/checks/A2_four_source_blocks.wl) | [A2_four_source_blocks.json](../../../data/kerr_point_particle_verification/results/A2_four_source_blocks.json) |
| A3：geodesic、两条 stress-tensor routes 与三个 Kinnersley components | [staged checks](../../../code/kerr_point_particle_verification/checks/a3/) | [A3 endpoint summary](../../../data/kerr_point_particle_verification/results/a3/a3_endpoint_summary.json)；`A3FullEndpointVerified=True` |
| A4：完整未投影点粒子方程到 fixed-$m$/$y$ endpoint | [staged checks](../../../code/kerr_point_particle_verification/checks/a4/) | [A4 endpoint summary](../../../data/kerr_point_particle_verification/results/a4/a4_endpoint_summary.json) |
| A5：common jets 与 LL/LJ/JJ/JL smooth-source interfaces | [staged checks](../../../code/kerr_point_particle_verification/checks/a5/) | [A5 endpoint summary](../../../data/kerr_point_particle_verification/results/a5/a5_endpoint_summary.json) |
| A6：一阶化、radial/angular characteristics、axis/projection contract | [targeted checks](../../../code/kerr_point_particle_verification/checks/a6/) | [A6 targeted summary](../../../data/kerr_point_particle_verification/results/a6/a6_targeted_summary.json)；FullEndpointClaimMade=false |

A1/A2 使用实数 $M,a,L,R,\theta$，$M,L,R>0$，$0<\theta<\pi$，并排除 $\Delta$ 及相关分母的零点。A3 另取非极端 Kerr、赤道 $Q=0$ inward geodesic、$\mathcal R\ge0$、$u^t,u^T>0$、$P+\sqrt{\mathcal R}\ne0$，并在 delta routes 中假设一个 local simple transverse root。A3 证明的是这些 assumptions 下的完整 local endpoint；不证明 global root uniqueness、global worldline monotonicity 或 turning-point differentiability。

证据强度并不完全相同。完整左端算符的 12-jet comparison、A2 source 的 end-to-end endpoint comparison、A3 invariant/BL 双路线到完整 projections endpoint 的 comparison，以及 A4 的完整 fixed-$m,y$ transformed-target comparison 是主要代数证据。A5 是 common-jets 与四块 generic jet-interface 的模块化闭包，不是单一完全展开生产源；A6 是 owner 要求范围内的 targeted symbolic checks，不是完整数值 endpoint。角算符同定义比较、rank-one projection factorization、Fourier finite matrix 和 ordering metadata 只是辅助一致性检查。

尚未关闭的物理边界为

$$
\boxed{
\text{度规号差与 Weyl scalar 的桥接},
\qquad
\text{有源方程的单位制约定},
\qquad
\text{Ripley Eq. (22) 的独立转录}
}.
$$

因此，本文当前最窄状态是：A2/A3 endpoints 已在各自 assumptions 下验证，A4 保留其既有条件代数结果；三者仍继承未关闭的 A1 convention/PDF boundary。

## V. 完整未分解 point-particle equation：冻结接口

本节冻结进入 fixed-$m$ 与 $y=-\cos\theta$ 之前的连续方程。坐标仍为 $(T,R,\theta,\Phi)$，field 仍为 Ripley regular-tetrad peeling variable $\psi_4(T,R,\theta,\Phi)$，没有作 Fourier projection、angular-coordinate change 或 Gaussian replacement。

### A. 完整左端

定义

$$
\boxed{
\mathcal O_{\rm R}^{(-2)}
\equiv\left.\mathcal O_{\rm R}\right|_{s=-2}
},
$$

其中 $\mathcal O_{\rm R}$ 是 Ripley Eq. (22) 在 $(T,R,\theta,\Phi)$ 中的完整零阶、一阶与二阶 operator；其逐项权威入口是 [附录 A 的完整 Eq. (22)](appendix/RIPLEY_OPERATOR.md#ripley-eq-22)。特别地，angular part 仍是未投影的 ${}_s\!\Delta_{S^2}$，且所有 $\partial_\Phi,\partial_T\partial_\Phi,\partial_R\partial_\Phi$ 均保留。

### B. Source operators 与三个 matter components

四块 source 使用

$$
\boxed{
\widehat{\mathcal L}_s
=\partial_\theta-i\csc\theta\,\partial_\Phi
-ia\sin\theta\,\partial_T+s\cot\theta
},
$$

$$
\boxed{
\widehat{\mathcal J}_+
=-\left(2+\frac{4M}{r}\right)\partial_T
-\frac{L^2}{r^2}\partial_R,
\qquad r=\frac{L^2}{R}
}.
$$

令

$$
\delta_p^{(3)}
=\delta(R-R_p(T))\delta(\theta-\theta_p)
\delta(\Phi-\Phi_p(T)),
$$

$$
\boxed{
\mathcal W_p^{(\theta)}
=\frac{\mu L^2}
{\Sigma_pr_p^2\sin\theta_p\,u_p^T}\delta_p^{(3)}
},
$$

以及 Kinnersley contractions

$$
\boxed{
\begin{aligned}
N_p={}&-\frac{\Delta_p}{2\Sigma_p}u_p^T
+\frac{\Delta_pr_p^2}{L^2\Sigma_p}
\left(1+\frac{2M}{r_p}\right)u_p^R
+\frac{a\Delta_p\sin^2\theta_p}{2\Sigma_p}u_p^\Phi,
\end{aligned}
}
$$

$$
\boxed{
\begin{aligned}
\bar M_p={}&\frac1{\sqrt2(r_p-ia\cos\theta_p)}
\Bigg[
ia\sin\theta_p\,u_p^T
-\frac{iar_p^2\sin\theta_p}{L^2}
\left(2+\frac{4M}{r_p}\right)u_p^R\\
&+\Sigma_pu_p^\theta
-i(r_p^2+a^2)\sin\theta_p\,u_p^\Phi
\Bigg].
\end{aligned}
}
$$

则三个完整 distributional components 是

$$
\boxed{
T_{nn}^{H,\rm pp}=\mathcal W_p^{(\theta)}N_p^2,
\qquad
T_{\bar mn}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_pN_p,
\qquad
T_{\bar m\bar m}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_p^2
}.
$$

对本项目的赤道 orbit，$\theta_p=\pi/2$、$u_p^\theta=0$、$\Sigma_p=r_p^2$，而 $u_p^T,u_p^R,u_p^\Phi,N_p,\bar M_p$ 由 Section III 的 horizon-regular geodesic formulas 给定。

### C. 四块 source 与完整方程

保持 outer-to-inner ordering，四块结构为

$$
\boxed{
\begin{aligned}
\widehat T_{H,1}^{\rm pp}
&=\widehat{\mathcal B}_{nn}^{(LL)}
[\mathcal W_p^{(\theta)}N_p^2],\\
\widehat T_{H,2}^{\rm pp}
&=\widehat{\mathcal B}_{\bar mn}^{(LJ)}
[\mathcal W_p^{(\theta)}\bar M_pN_p],\\
\widehat T_{H,3}^{\rm pp}
&=\widehat{\mathcal B}_{\bar m\bar m}^{(JJ)}
[\mathcal W_p^{(\theta)}\bar M_p^2],\\
\widehat T_{H,4}^{\rm pp}
&=\widehat{\mathcal B}_{\bar mn}^{(JL)}
[\mathcal W_p^{(\theta)}\bar M_pN_p].
\end{aligned}
}
$$

这里 $(LJ)$ 表示 outer $\widehat{\mathcal L}$ / inner $\widehat{\mathcal J}$，$(JL)$ 表示 outer $\widehat{\mathcal J}$ / inner $\widehat{\mathcal L}$；两块不能交换。四个 blocks 的全部 $\rho,\bar\rho,\Delta$ prefactors、operator nesting 和已经代入的 point-particle operands 逐项写在 [附录 C.6](appendix/POINT_PARTICLE_PROJECTIONS.md#c6-逐项代入后的-point-particle-four-blocks)，其 A2 functional 原型见 [附录 B](appendix/SASAKI_SOURCE.md#拉回后的-source)。

冻结的完整未分解方程是

$$
\boxed{
\mathcal O_{\rm R}^{(-2)}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}
\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}
}.
$$

依赖与求导范围必须按下列边界解释：

- **field point：** $\rho,\bar\rho,\Delta,\Sigma,r$ 与 source-operator coefficients 是 $(R,\theta)$ functions；方程外部的 $\Sigma/(\Delta^2R)$ 也是 field-point multiplier；
- **worldline：** $r_p,R_p,\theta_p,\Phi_p,u_p^A,N_p,\bar M_p$ 是 $T$ 的 functions，并在 particle support 上取值；
- **distribution：** $\mathcal W_p^{(\theta)}$ 包含相对于 $dR\,d\theta\,d\Phi$ 的三个 coordinate deltas；
- **operator action：** 每个 nested $\widehat{\mathcal L}_s$、$\widehat{\mathcal J}_+$ 都作用于其右侧方括号中的全部 background factors、worldline amplitudes 和 distributions。不得在求导前把 field-point factors 误当成 worldline constants，也不得交换 mixed-block ordering。

**冻结状态。** A2 已验证 arbitrary Kinnersley projection functions 从 BL source 到完整 hyperboloidal four-block endpoint；A3 已验证 geodesic、两条 stress-tensor routes 与三个具体 point-particle components，且 `A3FullEndpointVerified=True`。因此上式是进入下一步 fixed-$m/y$ transformation 的 verified unprojected algebraic interface。它仍继承 A1 的 signature/Weyl-scalar、sourced-unit 与 independent Ripley-PDF transcription boundaries；这里的“冻结”不表示 owner approval 或 full Gate A closure。

## VI. 固定 $m$ 与 $y=-\cos\theta$

本阶段严格从 Section V 冻结的完整未分解 equation 出发：

$$
\mathcal O_{\rm R}^{(-2)}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}
\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}.
$$

只依次作 azimuthal Fourier projection 与 $y=-\cos\theta$ coordinate transformation；不引入 Gaussian，也不展开 nested spatial derivatives。

### A. Fixed-$m$：$(T,R,\theta,\Phi)\to(T,R,\theta;m)$

采用 Fourier pair

$$
\boxed{
f=\sum_{m\in\mathbb Z}f_m(T,R,\theta)e^{im\Phi},
\qquad
f_m=\frac1{2\pi}\int_0^{2\pi}fe^{-im\Phi}\,d\Phi
},
$$

故

$$
\boxed{
\partial_\Phi\to im,
\qquad
\partial_\Phi^2\to-m^2
}.
$$

Point-particle periodic delta 同时投影：

$$
\boxed{
\delta_{2\pi}(\Phi-\Phi_p)
=\frac1{2\pi}\sum_m e^{im\Phi}e^{-im\Phi_p}
}.
$$

因此每个 source mode 含且只含一次 $(2\pi)^{-1}e^{-im\Phi_p(T)}$。Source angular operator 是

$$
\widehat{\mathcal L}_s^{[m,\theta]}
=\partial_\theta+m\csc\theta
-ia\sin\theta\,\partial_T+s\cot\theta,
$$

而 $\widehat{\mathcal J}_+$ 保持 Section V 的形式。定义

$$
\mathcal W_m^{(\theta)}
=\frac{\mu L^2}
{2\pi\Sigma_pr_p^2\sin\theta_p\,u_p^T}
e^{-im\Phi_p(T)}
\delta(R-R_p)\delta(\theta-\theta_p),
$$

$$
\boxed{
T_{nn,m}^{H,\theta}=\mathcal W_m^{(\theta)}N_p^2,
\quad
T_{\bar mn,m}^{H,\theta}=\mathcal W_m^{(\theta)}\bar M_pN_p,
\quad
T_{\bar m\bar m,m}^{H,\theta}=\mathcal W_m^{(\theta)}\bar M_p^2
}.
$$

把它们代入保持 A2 ordering 的四块后，得到完整中间 equation

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,\theta)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{(\theta),\rm dist}
}.
$$

$\mathcal O_{-2,m}^{(T,R,\theta)}$ 的全部 coefficients、fixed-$m$ angular operator 和四个 exact $\theta$-blocks 分别见 [附录 D.2.3](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md#a4-theta-operator)与 [附录 D.2.5](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md#a4-theta-blocks)。这一步只用 $\partial_\Phi\to im$ 和 periodic-delta coefficient，没有改变 $R,\theta$ distributions。

### B. Angular coordinate：$(T,R,\theta;m)\to(T,R,y;m)$

令

$$
\boxed{
y=-\cos\theta,
\qquad D=1-y^2,
\qquad
\sin\theta=\sqrt D,
\qquad
\partial_\theta=\sqrt D\,\partial_y
}.
$$

Measure、angular delta 和 tensor components 必须同时变换：

$$
\boxed{
\sqrt{-g_{TRy\Phi}}=\frac{\Sigma r^2}{L^2},
\qquad
\delta(\theta-\theta_p)=\sin\theta_p\,\delta(y-y_p),
\qquad y_p=-\cos\theta_p
},
$$

$$
u^y=\sin\theta\,u^\theta,
\qquad
v_y=\frac{v_\theta}{\sin\theta},
$$

$$
\boxed{
\Sigma=r^2+a^2y^2,
\qquad
\rho=(r+iay)^{-1},
\qquad
\bar\rho=(r-iay)^{-1}
}.
$$

因此 scalar contractions $N_p,\bar M_p$ 不变，而 $1/\sin\theta_p$ 与 angular-delta Jacobian 相消。最终公共 weight 与三个 components 是

$$
\boxed{
\mathcal W_m^{(y)}
=\frac{\mu L^2}
{2\pi\Sigma_pr_p^2u_p^T}
e^{-im\Phi_p(T)}
\delta(R-R_p(T))\delta(y-y_p(T))
},
$$

$$
\boxed{
T_{nn,m}^{H,y}=\mathcal W_m^{(y)}N_p^2,
\quad
T_{\bar mn,m}^{H,y}=\mathcal W_m^{(y)}\bar M_pN_p,
\quad
T_{\bar m\bar m,m}^{H,y}=\mathcal W_m^{(y)}\bar M_p^2
}.
$$

Angular operators 为

$$
\boxed{
\mathcal A_{s,m}[f]
=D f_{,yy}-2yf_{,y}
-\frac{(m-sy)^2}{D}f+sf
},
$$

$$
\boxed{
\widehat{\mathcal L}^{[m]}_s
=\sqrt D\,\partial_y
-ia\sqrt D\,\partial_T
+\frac{m-sy}{\sqrt D}
},
$$

而 $\widehat{\mathcal J}_+$ 不变。完整 field coefficients 见 [附录 D.4](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md#a4-y-operator)，四个保持 outer-to-inner ordering 的 exact blocks 见 [附录 D.5](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md#a4-y-source)。最终 pre-Gaussian equation 是

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,y)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{\rm dist}
}.
$$

### C. 分布积分意义、旋转轴正则性与证据状态

四块右端不是普通 pointwise functions，而是含 Dirac deltas 的 spacetime distributions。若

$$
F=C(T)\delta(R-R_p(T))\delta(y-y_p(T)),
$$

则对任意 compactly supported smooth test function $\chi(T,R,y)$，定义

$$
\langle F,\chi\rangle
=\int C(T)\chi(T,R_p(T),y_p(T))\,dT,
$$

$$
\boxed{
\left\langle\partial_R^j\partial_y^kF,\chi\right\rangle
=(-1)^{j+k}
\left\langle F,\partial_R^j\partial_y^k\chi\right\rangle
}.
$$

也就是说，nested $\partial_R,\partial_y$ 对 delta 的作用通过 integration by parts 定义，而不是在粒子位置作普通函数求导。每移动一个 derivative 都产生一个负号；nested operators 的 distributional transposes 必须反序作用。因此上述分布定义同时固定 delta-derivative signs 和两个 mixed blocks 的 operator ordering。这里保持这一积分定义，不展开巨大 delta-derivative expressions，也不消去 source 中的 $\partial_T$。

这里的“轴端”是 Kerr 旋转轴的两个 coordinate endpoints：

$$
y=-1\Longleftrightarrow\theta=0,
\qquad
y=+1\Longleftrightarrow\theta=\pi,
$$

不是复分析中的 poles。$\mathcal A_{s,m}$ 中 $(1-y^2)^{-1}$ 的表面奇异由 regular spin-weighted mode 在轴上的消失阶完成。对 $s=-2$，north/south 的 regular powers 是

$$
\boxed{
(p_-,p_+)_{m=2}=(0,2),
\qquad
(p_-,p_+)_{m=4}=(1,3)
}.
$$

即 regular mode 局部满足 $\psi_{4,m}\sim(1+y)^{p_-}(1-y)^{p_+}$ 乘以轴上 smooth function。详细 distributional definition 与 axis indicial derivation 见 [附录 D.6--D.7](appendix/W02_1_4_FIXED_M_Y_DERIVATION.md#a4-weak-poles)。这些 powers 与 Wigner-$d$ 检查是轴正则性的辅助证据；它们尚未独立闭合数值 angular boundary condition，因此不作为 A4 source transformation 的主结论。

**A4 证据状态。** Mathematica 从 Section V 的完整未投影周期点粒子方程 $E_0$ 出发，不调用理论的 $E_m^{(\theta)}$ 或 $E_m^{(y)}$ 作为生成 helper。Kernel 先执行 $\Phi$ 求导和 mode extraction，再执行 $y=-\cos\theta$、angular-delta Jacobian 与外层 $\Sigma$ 因子的变换。完整 $\theta$ 场算符、四个已代入 $\mathcal W_m^{(\theta)},N_p,\bar M_p$ 的源块和总方程 residual 均为 $0$；从同一 $E_0$ 直接到最终 $y$ 方程时，场算符、四个完整 $S_{m,A}^{\rm dist}$ 与端到端 residual 也均为 $0$。见 [A4a](../../../data/kerr_point_particle_verification/results/a4/a4a_fixed_m_theta.json)、[A4b](../../../data/kerr_point_particle_verification/results/a4/a4b_y_field_source.json)和 [汇总](../../../data/kerr_point_particle_verification/results/a4/a4_endpoint_summary.json)。因此 `A4FullEndpointVerified=True`。该结论使用 $R>0$、$-1<y<1$、$\Delta\ne0$ 且粒子不位于坐标轴端点；轴上连续正则性幂次已作辅助检查，但数值 angular closure 仍开放。A4 继承 A1--A3 的 assumptions，并排除 Gaussian、off-worldline extension、time jets 与 discretization。

## VII. Gaussian、radial extension 与 time jets

A4 的冻结输入是四个 exact blocks

$$
S_{m,A}^{\rm dist},
\qquad A=1,\ldots,4,
$$

及其完整 equation

$$
\mathcal O_{-2,m}^{(T,R,y)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{\rm dist}.
$$

本阶段只把 exact coordinate deltas Gaussianized，并解析消去四块中的全部 coordinate $\partial_T$；$\partial_R,\partial_y$ 保持原 ordering，留给有限差分 matrices。

### A. Gaussian 与 radial off-worldline extension

取固定、正的 widths $\sigma_R,\sigma_y$，相对于 $dR\,dy$ 定义 full-line normalized kernels

$$
G_R=\frac{e^{-(R-R_p)^2/(2\sigma_R^2)}}
{\sqrt{2\pi}\sigma_R},
\qquad
G_y=\frac{e^{-(y-y_p)^2/(2\sigma_y^2)}}
{\sqrt{2\pi}\sigma_y},
$$

$$
\boxed{
\delta(R-R_p)\delta(y-y_p)\longrightarrow G_RG_y
}.
$$

A4 已处理 angular Jacobian，此处不再加入 $\sin\theta$。定义

$$
\widehat\Delta=R^2\Delta,
\qquad
\widehat\Sigma=R^2\Sigma,
\qquad
\widehat N_p=\frac{N_p}{\widehat\Delta_p}
=-\frac{K_p}{2L^4},
$$

其中

$$
K=\frac{P-\mathcal V}{\Delta}
=\frac{C}{P+\mathcal V},
\qquad
\mathcal V=\sqrt{\mathcal R}.
$$

采用 compactified-$\Delta$ radial extension

$$
\boxed{
F_{nn}^{\rm ext}
=\widehat\Delta(R)^2\widetilde Q_{nn}G_RG_y,
\quad
F_{\bar mn}^{\rm ext}
=\widehat\Delta(R)\widetilde Q_{\bar mn}G_RG_y,
\quad
F_{\bar m\bar m}^{\rm ext}
=Q_{\bar m\bar m}G_RG_y
}.
$$

在 $R=R_p$ 上它恢复 exact projection coefficients；在 Gaussian replacement 前与 A4 source 是同一 distribution。有限 widths 下它是一项明确但非唯一的 regularization choice。完整 field-point/worldline factor boundary 见 [附录 E.6](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-extension)。

### B. 已闭合的二阶 time jets

由 regular geodesic 解析生成

$$
\left\{
r_p,R_p,\Phi_p,u_p^T,u_p^R,u_p^\Phi,
\widehat N_p,\bar M_p
\right\}^{(k)},
\qquad k=0,1,2.
$$

例如

$$
\dot R_p=\frac{u_p^R}{u_p^T},
\qquad
\ddot R_p
=\dot r_p\frac{u^R_{p,r}u_p^T-u_p^Ru^T_{p,r}}
{(u_p^T)^2},
$$

$$
\dot\Phi_p=\frac{u_p^\Phi}{u_p^T},
\qquad
\ddot\Phi_p
=\dot r_p\frac{u^\Phi_{p,r}u_p^T-u_p^\Phi u^T_{p,r}}
{(u_p^T)^2}.
$$

全部 radial derivatives、velocity jets、$\widehat N_p,\bar M_p$ jets 与三个实际 extension amplitudes 的显式二阶 formulas 见 [附录 E.4--E.5](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-worldline-jets)。没有把 dot quantities 留作外部输入。

对 exact 或 tilde amplitude，

$$
Q_{ab}=A_{ab}e^{-im\Phi_p},
$$

$$
\dot Q_{ab}
=e^{-im\Phi_p}(\dot A_{ab}-im\dot\Phi_pA_{ab}),
$$

$$
\ddot Q_{ab}
=e^{-im\Phi_p}
\left[
\ddot A_{ab}-2im\dot\Phi_p\dot A_{ab}
-im\ddot\Phi_pA_{ab}-m^2\dot\Phi_p^2A_{ab}
\right].
$$

对 $Z(T)G_RG_y$，固定 widths 给出

$$
\boxed{
f_Z^{(0)}=ZG,
\quad
f_Z^{(1)}=(\dot Z+Z\Gamma_1)G,
\quad
f_Z^{(2)}=(\ddot Z+2\dot Z\Gamma_1+Z\Gamma_2)G
},
$$

其中赤道 $y_p=0$ 时

$$
\Gamma_1=\frac{(R-R_p)\dot R_p}{\sigma_R^2},
$$

$$
\Gamma_2
=\frac{(R-R_p)^2\dot R_p^2}{\sigma_R^4}
+\frac{(R-R_p)\ddot R_p-\dot R_p^2}{\sigma_R^2}.
$$

完整 projection-jet triples 见 [附录 E.7](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-gaussian-jets)。

### C. 四个数值源块

对已知 jet sequence $X^{(0)},X^{(1)},X^{(2)}$，定义只含 spatial derivatives 的 operators

$$
c=\sqrt{1-y^2},
\qquad
\ell_s=\frac{m-sy}{c},
\qquad
j_T=-\left(2+\frac{4MR}{L^2}\right),
\qquad
j_R=-\frac{R^2}{L^2},
$$

$$
\boxed{
\mathfrak L_s[X]^{(k)}
=c\,\partial_yX^{(k)}-iac\,X^{(k+1)}
+\ell_sX^{(k)},
\qquad k=0,1
},
$$

$$
\boxed{
\mathfrak J[X]^{(k)}
=j_TX^{(k+1)}+j_R\partial_RX^{(k)},
\qquad k=0,1
}.
$$

它们是把 $\partial_TX^{(k)}=X^{(k+1)}$ 解析代入 $\widehat{\mathcal L}_s^{[m]}$ 与 $\widehat{\mathcal J}_+$ 后的记账形式。逐块得到

$$
\boxed{
S_{m,A}^{\rm num}
=-16\pi R\widehat\Sigma(R,y)\,\mathfrak H_A,
\qquad A=1,2,3,4
},
$$

其中

$$
\begin{array}{c|c|c}
A & \text{matter operand} & \text{outer-to-inner ordering}\\ \hline
1 & \widetilde f_{nn}^{(0,1,2)} & LL\\
2 & \widetilde f_{\bar mn}^{(0,1,2)} & LJ\\
3 & f_{\bar m\bar m}^{(0,1,2)} & JJ\\
4 & \widetilde f_{\bar mn}^{(0,1,2)} & JL
\end{array}
$$

且

$$
\boxed{
S_m^{\rm num}=\sum_{A=1}^4S_{m,A}^{\rm num}
}.
$$

四个 $\mathfrak H_A$ 的完整 nested definitions 见 [附录 E.8](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-four-blocks)。它们不含 Dirac delta、不含 coordinate $\partial_T$、不含未定义 jets；只保留 $\partial_R,\partial_y$ 作用于已知 smooth operands。在 $0<R<R_H$、$|y|<1$ 内，compactified-$\Delta$ regrouping 还使四块分别无 horizon/SCRI+ radial Laurent 负幂。

### D. A5 状态与开放边界

A5 的 exact scalar inputs、projection-jet outputs 和四个逐块 residual targets 见 [附录 E.10](appendix/W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-contract)。MMA 必须从各自 A4 block 直接执行 $\partial_T$，再与对应 $S_{m,A}^{\rm num}$ 比较；四块全部为零后才比较 total source。

**A5 验证状态。** Mathematica 从共同 scalar definitions 直接得到 worldline、$\widehat N_p,\bar M_p$、三个 extension amplitudes、phase 与 moving-Gaussian 的 $0,1,2$ 阶 jets；全部 residual 为 $0$。$G_R,G_y$ 的 full-line normalization 和 centered first moments residual 为 $0$，compactified-$\Delta$ extension 在 worldline 恢复 A4 coefficients 的四项 residual 也为 $0$。随后从四个 A4 blocks 分别直接执行全部 $\partial_T$，与 $S_{m,A}^{\rm num}$ 比较得到

$$
\operatorname{Res}(S_{m,1}^{\rm num},S_{m,2}^{\rm num},
S_{m,3}^{\rm num},S_{m,4}^{\rm num})=(0,0,0,0),
\qquad
\operatorname{Res}S_m^{\rm num}=0.
$$

生成路线不调用理论的 $\mathfrak L_s,\mathfrak J,\mathfrak H_A$ helpers，且输出不含 Dirac delta、未求值 $\partial_T$ 或未定义 jets；LL、LJ、JJ、JL 空间次序保持不变。见 [A5 common jets](../../../data/kerr_point_particle_verification/results/a5/a5a_common_time_jets.json)、[四块结果](../../../data/kerr_point_particle_verification/results/a5/)和 [汇总](../../../data/kerr_point_particle_verification/results/a5/a5_endpoint_summary.json)。因此 `A5FullEndpointVerifiedWithinInteriorScope=True`。Full-line Gaussian 在有限 computational domain 上不严格单位归一；普通 $G_y$ 的 axis tails 尚无 projection-specific spin-frame completion；$\mathcal V=0$ turning point 需要另一参数化；有限-width convergence 属于 W02.3。这些边界没有被 time-algebra residual 关闭。

供物理复核和后续 Python finite-difference 实现直接使用的四块逐层公式，集中列在 [A5 四块数值源公式页](appendix/a5_four_block_source_terms.md)。

## VIII. A6 final PDE 与数值接口

本节把已验证的 interior source $S_m^{\rm num}$ 接到完整 Ripley operator，并给出 W02.3 所需的连续接口。详细 coefficients、axis derivation、stage data flow 和工程边界集中在 [A6 数值接口附录](appendix/a6_numerical_interface.md)。

### A. Post-ISCO inputs 与唯一轨迹

由 Kerr prograde ISCO 得到 $(r_{\rm ISCO},E_{\rm ISCO},L_{z,\rm ISCO})$ 后，实际 plunge input 写为

$$
\boxed{
r_0=r_{\rm ISCO}-\delta r,
\qquad
E_{\rm p}=E_{\rm ISCO}-\delta E,
\qquad
L_{z,\rm p}=L_{z,\rm ISCO}-\delta L_z
},
$$

并要求

$$
r_+<r_0<r_{\rm ISCO},
\qquad
\mathcal R(r_0)\ge0,
\qquad
E_{\rm p}-\Omega_HL_{z,\rm p}>0.
$$

Transition prescription 决定三个 deltas，但不是上述 inequalities 的代数结果。数值只积分一条轨迹

$$
\boxed{
z_p(T)=\bigl(R_p(T),0,\Phi_p(T)\bigr),
\qquad
\dot R_p=\frac{u_p^R}{u_p^T},
\qquad
\dot\Phi_p=\frac{u_p^\Phi}{u_p^T}
},
$$

初值为 $R_p(0)=L^2/r_0$、$\Phi_p(0)=0$。$r_p=L^2/R_p$ 只用于评价 BL radial functions，不形成第二条 trajectory。

### B. 二阶 PDE 与连续 P/Q 系统

记 $\psi=\psi_{4,m}$，完整 sourced equation regroup 为

$$
\boxed{
A\psi_{,TT}+B\psi_{,TR}+C\psi_{,RR}
-\mathcal A_{-2,m}[\psi]
+D\psi_{,T}+E_R\psi_{,R}+F\psi
=S_m^{\rm num}
}.
$$

$A,B,C,D,E_R,F$ 的完整 expressions 由附录 D.4 取 $s=-2$ 得到，并逐项列在 [A6 附录](appendix/a6_numerical_interface.md#a6-pde)。按 `P13` Eq. (23) 定义

$$
\boxed{
P=A\psi_{,T}+B\psi_{,R}+D\psi,
\qquad
Q=\psi_{,R}
},
$$

则

$$
\boxed{
\psi_{,T}=\frac{P-BQ-D\psi}{A}
},
$$

$$
\boxed{
P_{,T}
=S_m^{\rm num}-CQ_{,R}
+\mathcal A_{-2,m}[\psi]-E_RQ-F\psi
},
$$

$$
\boxed{
Q_{,T}=\partial_R
\left(\frac{P-BQ-D\psi}{A}\right)
}.
$$

连续 constraint 是 $Q-\psi_{,R}=0$。Method-of-lines baseline 只演化 $(P,\psi)$，每个 RK stage 重算

$$
Q=D_R^{(1)}\psi,
\qquad
Q_{,R}=D_R^{(2)}\psi.
$$

初始 $(\psi_0,\Pi_0)$ 转为

$$
\boxed{
P_0=A\Pi_0+B D_R^{(1)}\psi_0+D\psi_0
}.
$$

每个 RK4 stage 都重新计算 radial/angular derivatives、采样 $z_p(T_s)$、生成 A5 jets 和四个 $S_{m,A}^{\rm num}(T_s)$，再评价 $(\dot P_s,\dot\psi_s)$；不能只在完整 step 起点更新 source。

### C. Characteristics 与 radial directions

Radial principal speeds 是

$$
\boxed{
v_R^\pm
=\frac{B/2\pm\sqrt{B^2/4-AC}}{A}
}.
$$

在 $A>0$ 的 nonextremal accepted domain，

$$
\boxed{
(v_R^+,v_R^-)_{R=0}
=\left(0,-\frac{2L^2}{A(0,y)}\right)
},
$$

$$
\boxed{
(v_R^+,v_R^-)_{R=R_H}
=\left(\frac{4MR_H}{A(R_H,y)},0\right)
}.
$$

SCRI+ 与 future horizon 均只有 tangent/outgoing characteristics，没有指向 exterior computational domain 的 incoming physical field data。One-sided rows 是数值 closure，不是额外 boundary data。

### D. Axis regularity 与局部 closure

对 $s=-2$，axis indicial equations 给出

$$
\boxed{
p_-=\frac{|m-2|}{2},
\qquad
p_+=\frac{|m+2|}{2}
},
$$

其中 $y=-1$ 是 north axis、$y=+1$ 是 south axis。Baseline 采用 open angular grid，不作 global field redefinition；局部写

$$
\psi=f_\pm g_\pm,
\qquad
f_-=(1+y)^{p_-},
\qquad
f_+=(1-y)^{p_+},
$$

并用 factor-aware weights

$$
\boxed{
\sum_jw_{ij}^{(n,\pm)}f_\pm(y_j)(y_j-y_i)^k
=\left.
\partial_y^n\left[f_\pm(y)(y-y_i)^k\right]
\right|_{y_i},
\quad k=0,\ldots,N_{\rm st}-1
}.
$$

$p_\pm>0$ 时 endpoint field 为零；$p_\pm=0$ 时由 smooth $g_\pm$ local polynomial 重建。完整 angular operator 必须先作用于 $f_\pm g_\pm$ 再取 finite axis limit，不能分别评价 $(1-y^2)^{-1}$ terms。具体 nodes、weights 和 finite-limit rows 等待 W02.3 验证。

### E. SCRI+ projection 与 Gaussian strategy

在 $R=0$、$T=u$，

$$
\boxed{
\psi_{4,\ell m}(T)
=2\pi\int_{-1}^{1}
{}_{-2}\!\bar Y_{\ell m}(\arccos(-y),0)
\psi_{4,m}(T,0,y)\,dy
}.
$$

离散 projection matrix 是

$$
\boxed{
\mathsf P_j^{(\ell m)}
=2\pi w_j{}_{-2}\!\bar Y_{\ell m}
(\arccos(-y_j),0),
\qquad
\psi_{4,\ell m}\approx\sum_j\mathsf P_j^{(\ell m)}\psi_j
}.
$$

Gaussian production baseline 仅作为工程策略取

$$
\boxed{
\sigma_R=4\Delta R,
\qquad
\sigma_y=4\Delta y
}.
$$

它定义 coupled grid--Gaussian refinement path，不是连续理论唯一宽度。

**A6 targeted-check 状态。** Mathematica 已完成三组必要检查：[一阶化](../../../data/kerr_point_particle_verification/results/a6/a6_first_order.json) 的 A--F、P/Q 与 source placement residual 全为 $0$；[characteristics](../../../data/kerr_point_particle_verification/results/a6/a6_characteristics.json) 的 radial/angular principal roots、SCRI+/horizon limits 与 directions 通过；[axis/projection](../../../data/kerr_point_particle_verification/results/a6/a6_axis_projection.json) 的 indicial/Laurent、toy factor-aware basis 与 analytic projection residual 通过。[汇总](../../../data/kerr_point_particle_verification/results/a6/a6_targeted_summary.json) 给出 `AllTargetedChecksPassed=True` 且 `FullEndpointClaimMade=False`。Owner 接受当前连续 contract 和 targeted symbolic scope，但明确不把 boundary-adjacent finite-difference rows 视为已理解或通过。Transition、startup、production stencils、CFL/stability、axis finite-limit rows、quadrature choice、output cadence、coupled convergence 和 production admission 仍由 W02.3 关闭。

## IX. 当前边界

- A1 signature/Weyl-scalar、sourced-unit convention 与 Ripley PDF 独立转录仍开放。
- A2/A3 endpoints 与 A4 full endpoint 已在记录的 assumptions 下验证，并继承 A1 boundary。
- A5 interior time-algebra endpoint 已验证；finite-domain normalization、angular-axis completion、turning branch 和 width convergence 分开保持开放。
- A6 continuous contract 与三组 targeted checks 已获 owner 接受且未作 full-endpoint claim；boundary finite-difference closure、工程稳定性和收敛性仍开放。
- Full Gate A 尚未关闭，W02.3 数值实现不能把本文件未验证部分当作冻结接口。
