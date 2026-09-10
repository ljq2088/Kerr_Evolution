# 附录 E：Gaussian、off-worldline extension 与解析时间 jets

- 任务：T006，W02.1.5
- 状态：A5 interior time-algebra endpoint 已验证；A4 full endpoint 已验证并继承 A1；finite-domain/axis/turning/width 边界开放
- 条件输入：[W02.1.4 fixed-$m$/$y$ 候选](W02_1_4_FIXED_M_Y_DERIVATION.md)
- 上游：[DERIVATION.md](../DERIVATION.md)及相关 theory appendices
- 范围：连续理论；不含 Mathematica 执行、有限差分、grid、CFL、RK4 或数值宽度选择
- 证据边界：本文公式由 A4 exact blocks 与赤道 Kerr geodesic 显式推导；目标仓库 C01 的当前 source/time-jet 文件只用于定位 compactified-$\Delta$ 候选和交叉阅读，历史代码及其历史 residual 不能替代理论推导或本仓库 A5 验证

## E.1 逻辑链与当前结论

```text
A4 exact distributional four-block equation（条件输入）
  -> 选择相对于 dR dy 的 Gaussian 近似恒等核
  -> 明确 worldline / field-point factors
  -> 用 compactified-Delta 指定有限宽 radial extension
  -> 在 distribution level 证明与 A4 同一
  -> 用赤道 Kerr geodesic 构造 R_p、Phi_p 与 contraction jets
  -> 对 phase、amplitudes、moving Gaussian 求 0/1/2 阶 T jets
  -> 逐块解析执行所有 nested partial_T
  -> 每块只留下作用在光滑 interior operands 上的 partial_R、partial_y
  -> 分别检查 horizon 与 SCRI+ 的 Laurent 阶
```

在 A4 及其 A1--A3 前提成立时，本文得到一套明确的 **radial compactified-$\Delta$ conditional extension**。它在 $R=R_p$ 恢复 exact coefficient，并允许四块 source 分别写成不含 $\widehat\Delta^{-1}$ 的形式；由解析幂次计数，每块在 future horizon 与 SCRI+ 均无 Laurent 负幂，不依赖块间相消。

但还有一个会影响闭区间数值接口的开放问题：普通 $G_{\sigma_y}(y-y_p)$ 在 $[-1,1]$ 上不严格单位归一，并且其非零 axis tails 不自动提供 source intermediate operands 所需的 spin-frame axis completion。本文因此只把最终四块称为 $|y|<1$ 上的条件光滑源；在 projection-specific angular extension 或不含 axis endpoints 的 representation 被数值接口固定前，不声称已经得到闭区间 $y\in[-1,1]$ 上的最终数值源。A4 的 continuum field axis powers 是必要的下游边界信息，但它们本身不足以唯一决定三个 matter projections 的 off-worldline angular extension。

## E.2 A4 输入与记号

A4 full endpoint 已在继承 A1--A3 输入和记录的 assumptions 下通过 staged Mathematica closure；见 [A4 summary](../../../../data/kerr_point_particle_verification/results/a4/a4_endpoint_summary.json)。本文把其 boxed equation

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,y)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{\rm dist}
}
$$

作为冻结输入，不重新证明 Fourier phase、$1/(2\pi)$、angular Jacobian、四块顺序或整体 source sign。

使用

$$
r=\frac{L^2}{R},\qquad
\Delta=r^2-2Mr+a^2,\qquad
\Sigma=r^2+a^2y^2,
$$

并定义在 compactified coordinate 中为多项式的量

$$
\boxed{
\widehat\Delta(R)=R^2\Delta
=L^4-2ML^2R+a^2R^2
},
$$

$$
\boxed{
\widehat\Sigma(R,y)=R^2\Sigma
=L^4+a^2R^2y^2
}.
$$

于是 A4 的外部 field-point prefactor 可写为

$$
-\frac{16\pi\Sigma}{\Delta^2R}
=-\frac{16\pi R\widehat\Sigma}{\widehat\Delta^2}.
$$

另记

$$
c(y)=\sqrt{1-y^2},\qquad
\ell_s(y)=\frac{m-sy}{c(y)},
$$

$$
j_T(R)=-\left(2+\frac{4MR}{L^2}\right),
\qquad
j_R(R)=-\frac{R^2}{L^2},
$$

$$
\rho=\frac{R}{L^2+iaRy},\qquad
\bar\rho=\frac{R}{L^2-iaRy}.
$$

所有带下标 $p$ 的量均在 worldline 上求值；不带 $p$ 的 $R,y$ 是 field point。

## E.3 Gaussian replacement、固定 widths 与 measure

对固定但尚未数值选定的 $\sigma_R,\sigma_y>0$，定义 full-line normalized kernels

$$
G_R(T,R)=\frac1{\sqrt{2\pi}\sigma_R}
\exp\left[-\frac{(R-R_p(T))^2}{2\sigma_R^2}\right],
$$

$$
G_y(T,y)=\frac1{\sqrt{2\pi}\sigma_y}
\exp\left[-\frac{(y-y_p(T))^2}{2\sigma_y^2}\right].
$$

它们按 coordinate measure 归一：

$$
\int_{-\infty}^{\infty}G_R\,dR=1,
\qquad
\int_{-\infty}^{\infty}G_y\,dy=1.
$$

因此条件替换是

$$
\boxed{
\delta(R-R_p)\delta(y-y_p)
\longrightarrow G_R G_y
}.
$$

这里没有 $\sin\theta$、$1/\sin\theta$ 或其他 angular Jacobian。A4 已经使用

$$
\sqrt{-g_{TRy\Phi}}=\frac{\Sigma r^2}{L^2},
\qquad
\delta(\theta-\theta_p)=\sin\theta_p\delta(y-y_p)
$$

消去了原来的 angular sine；在本阶段再次加入会重复计算 Jacobian。

这一定义是对开放域内 delta sequence 的连续理论选择。对实际有限域 $R\in[0,R_H]$、$y\in[-1,1]$，普通 Gaussian 的截断积分一般不严格等于 $1$。是否采用 full-line kernel、domain-renormalized kernel，或把源延拓到 horizon 内侧，是 W02.3 必须固定的数值接口；本文不静默选择 domain renormalization，因为其随 $R_p$ 变化的 normalization 会产生额外 time jets，而当前 A5 time-algebra target 不包含这些项。

<a id="a5-worldline-jets"></a>

## E.4 赤道 Kerr worldline 与完整 radial jets

### E.4.1 BL geodesic 与 horizon-finite velocity

令

$$
B=L_z-aE,
\qquad
P(r)=E(r^2+a^2)-aL_z,
$$

$$
C(r)=r^2+B^2,
\qquad
\mathcal R(r)=P(r)^2-\Delta(r)C(r),
\qquad
\mathcal V(r)=\sqrt{\mathcal R(r)}.
$$

本文采用 future-directed inward branch

$$
u^r=-\frac{\mathcal V}{r^2}.
$$

BL 其余分量为

$$
u^t=\frac1{r^2}
\left[aB+\frac{(r^2+a^2)P}{\Delta}\right],
$$

$$
u^\phi=\frac1{r^2}
\left[B+\frac{aP}{\Delta}\right].
$$

直接在 horizon 附近使用 $P/\Delta$ 会产生数值上不必要的 divergent pieces。由

$$
(P-\mathcal V)(P+\mathcal V)
=P^2-\mathcal V^2=\Delta C
$$

定义

$$
\boxed{
K(r)=\frac{P-\mathcal V}{\Delta}
=\frac{C}{P+\mathcal V}
}.
$$

只要 future-horizon branch 满足 $P(r_+)+\mathcal V(r_+)\ne0$，$K$ 有有限极限。把 BL velocity 代入坐标 Jacobian 后得到

$$
\boxed{
u^T=\frac1{r^2}
\left[
aB+(r^2+a^2)K
+\left(2+\frac{4M}{r}\right)\mathcal V
\right]
},
$$

$$
\boxed{
u^R=\frac{L^2\mathcal V}{r^4},
\qquad
u^\Phi=\frac{B+aK}{r^2}
}.
$$

这三个量在上述 branch assumptions 下 horizon-finite。若 $\mathcal V=0$（精确 radial turning point）或 $P+\mathcal V=0$，下述 radial jets 分支失效；本文不跨过该点作延拓。

### E.4.2 显式 radial jets

所有 worldline time jets 可由关于 $r$ 的解析导数产生。先写

$$
P'=2Er,\quad P''=2E,
\qquad
\Delta'=2(r-M),\quad\Delta''=2,
$$

$$
C'=2r,\qquad C''=2,
$$

$$
\mathcal R'
=2PP'-\Delta'C-\Delta C',
$$

$$
\mathcal R''
=2(P'^2+PP'')-\Delta''C-2\Delta'C'-\Delta C''.
$$

在 $\mathcal V>0$ 的 branch 上，

$$
\mathcal V'=\frac{\mathcal R'}{2\mathcal V},
\qquad
\mathcal V''=\frac{\mathcal R''}{2\mathcal V}
-\frac{(\mathcal R')^2}{4\mathcal V^3}.
$$

令 $H=P+\mathcal V$，则

$$
K=\frac{C}{H},
$$

$$
K'=\frac{C'}H-\frac{CH'}{H^2},
$$

$$
K''=\frac{C''}H
-\frac{2C'H'}{H^2}
-\frac{CH''}{H^2}
+\frac{2C(H')^2}{H^3},
$$

其中 $H'=P'+\mathcal V'$、$H''=P''+\mathcal V''$。

定义

$$
\eta(r)=2+\frac{4M}{r},
\qquad
\eta'=-\frac{4M}{r^2},
\qquad
\eta''=\frac{8M}{r^3},
$$

$$
F_T=aB+(r^2+a^2)K+\eta\mathcal V,
\qquad
F_\Phi=B+aK.
$$

所需 derivatives 全部为

$$
F_T'
=2rK+(r^2+a^2)K'
+\eta'\mathcal V+\eta\mathcal V',
$$

$$
F_T''
=2K+4rK'+(r^2+a^2)K''
+\eta''\mathcal V+2\eta'\mathcal V'
+\eta\mathcal V'',
$$

$$
F_\Phi'=aK',
\qquad
F_\Phi''=aK''.
$$

若 $u=F/r^2$，则

$$
u_{,r}=\frac{F'}{r^2}-\frac{2F}{r^3},
\qquad
u_{,rr}=\frac{F''}{r^2}-\frac{4F'}{r^3}
+\frac{6F}{r^4}.
$$

分别对 $(F_T,F_\Phi)$ 使用上式，得到 $u^T_{,r},u^T_{,rr}$ 与 $u^\Phi_{,r},u^\Phi_{,rr}$。径向分量显式为

$$
u^R_{,r}=L^2\left(
\frac{\mathcal V'}{r^4}-\frac{4\mathcal V}{r^5}
\right),
$$

$$
u^R_{,rr}=L^2\left(
\frac{\mathcal V''}{r^4}
-\frac{8\mathcal V'}{r^5}
+\frac{20\mathcal V}{r^6}
\right).
$$

### E.4.3 $R_p$ 与 $\Phi_p$ 的 time jets

点表示 $d/dT$。由 $dT/d\tau=u^T$，

$$
\boxed{
\dot r_p=\frac{u^r_p}{u^T_p}
=-\frac{r_p^2u^R_p}{L^2u^T_p}
}.
$$

令

$$
v(r)\equiv-\frac{r^2u^R}{L^2u^T}=\frac{u^r}{u^T}.
$$

其 radial derivative 是

$$
v'(r)
=-\frac1{L^2}
\left[
\frac{2ru^R+r^2u^R_{,r}}{u^T}
-\frac{r^2u^Ru^T_{,r}}{(u^T)^2}
\right].
$$

因此

$$
\ddot r_p
=v_pv'_{p}.
$$

更直接地，

$$
\boxed{
R_p=\frac{L^2}{r_p},
\qquad
\dot R_p=\frac{u^R_p}{u^T_p}
},
$$

$$
\boxed{
\ddot R_p
=\dot r_p
\frac{u^R_{p,r}u^T_p-u^R_pu^T_{p,r}}
(u^T_p)^2}
}.
$$

相位轨道满足

$$
\boxed{
\dot\Phi_p=\frac{u^\Phi_p}{u^T_p}
},
$$

$$
\boxed{
\ddot\Phi_p
=\dot r_p
\frac{u^\Phi_{p,r}u^T_p-u^\Phi_pu^T_{p,r}}
(u^T_p)^2}
}.
$$

任意只依赖 $r_p$ 的量 $X$ 均有

$$
\dot X=X_{,r}\dot r_p,
\qquad
\ddot X=X_{,rr}\dot r_p^2+X_{,r}\ddot r_p.
$$

特别地，对 $A\in\{T,R,\Phi\}$，

$$
\boxed{
\dot u_p^A=u^A_{p,r}\dot r_p,
\qquad
\ddot u_p^A=u^A_{p,rr}\dot r_p^2
+u^A_{p,r}\ddot r_p
}.
$$

至此 $r_p,R_p,\Phi_p$ 和三个 regular velocity components 的 $0,1,2$ 阶 $T$ jets 都只依赖 $(M,a,L,E,L_z,r_p)$ 与本节显式 radial functions，没有未定义 dot quantities。

<a id="a5-amplitude-jets"></a>

## E.5 Contractions、amplitudes 与 phase jets

### E.5.1 Exact contractions 与 regularized radial contraction

A3 的 $N_p,\bar M_p$ 赤道公式给出

$$
N_p=-\frac{P_p-\mathcal V_p}{2r_p^2}
=-\frac{\Delta_pK_p}{2r_p^2},
$$

$$
\bar M_p=\frac{i(aE-L_z)}{\sqrt2r_p}
=-\frac{iB}{\sqrt2r_p}.
$$

由于

$$
\widehat\Delta_p=R_p^2\Delta_p
=\frac{L^4\Delta_p}{r_p^2},
$$

定义

$$
\boxed{
\widehat N_p=\frac{N_p}{\widehat\Delta_p}
=-\frac{K_p}{2L^4}
}.
$$

这个等式给出了 horizon 上 $0/0$ 的唯一连续值，不需要数值相除。它也可由 regular hyperboloidal contraction 直接写成

$$
\widehat N_p
=-\frac{u^T_p}{2L^4}
+\frac{r_p^2+2Mr_p}{L^6}u^R_p
+\frac{a}{2L^4}u^\Phi_p.
$$

两条表达式在 A3 条件输入成立时等价。

径向 derivatives 为

$$
\widehat N_{,r}=-\frac{K'}{2L^4},
\qquad
\widehat N_{,rr}=-\frac{K''}{2L^4},
$$

$$
\bar M_{,r}=\frac{iB}{\sqrt2r^2},
\qquad
\bar M_{,rr}=-\frac{2iB}{\sqrt2r^3}.
$$

所以完整 time jets 是

$$
\boxed{
\dot{\widehat N}_p
=-\frac{K'_p}{2L^4}\dot r_p,
\qquad
\ddot{\widehat N}_p
=-\frac1{2L^4}
\left(K''_p\dot r_p^2+K'_p\ddot r_p\right)
},
$$

$$
\boxed{
\dot{\bar M}_p
=\frac{iB}{\sqrt2r_p^2}\dot r_p,
\qquad
\ddot{\bar M}_p
=-\frac{2iB}{\sqrt2r_p^3}\dot r_p^2
+\frac{iB}{\sqrt2r_p^2}\ddot r_p
}.
$$

为同时追踪 exact amplitudes，记 $h_p=\widehat\Delta(R_p)$。则

$$
N_p=h_p\widehat N_p,
$$

$$
\dot h_p=\widehat\Delta_{,R}(R_p)\dot R_p,
$$

$$
\ddot h_p
=\widehat\Delta_{,RR}(R_p)\dot R_p^2
+\widehat\Delta_{,R}(R_p)\ddot R_p,
$$

其中

$$
\widehat\Delta_{,R}=-2ML^2+2a^2R,
\qquad
\widehat\Delta_{,RR}=2a^2.
$$

因此

$$
\dot N=\dot h_p\widehat N+h_p\dot{\widehat N},
$$

$$
\ddot N
=\ddot h_p\widehat N+2\dot h_p\dot{\widehat N}
+h_p\ddot{\widehat N}.
$$

### E.5.2 Common weight 与完整 amplitude jets

赤道 $y_p=0$ 时，A4 projection 的 common worldline weight 是

$$
\boxed{
A_0=\frac{\mu L^2}{2\pi r_p^4u^T_p}
}.
$$

定义

$$
\Lambda_0=4\frac{\dot r_p}{r_p}
+\frac{\dot u^T_p}{u^T_p},
$$

$$
\dot\Lambda_0
=4\left(\frac{\ddot r_p}{r_p}
-\frac{\dot r_p^2}{r_p^2}\right)
+\frac{\ddot u^T_p}{u^T_p}
-\left(\frac{\dot u^T_p}{u^T_p}\right)^2.
$$

则

$$
\dot A_0=-\Lambda_0A_0,
\qquad
\ddot A_0=(\Lambda_0^2-\dot\Lambda_0)A_0.
$$

这里

$$
\dot u^T_p=u^T_{p,r}\dot r_p,
\qquad
\ddot u^T_p=u^T_{p,rr}\dot r_p^2
+u^T_{p,r}\ddot r_p.
$$

对于任意 $X,Y$，令 $A=A_0XY$。其完整 jets 是

$$
A^{(0)}=A_0XY,
$$

$$
\dot A=\dot A_0XY+A_0(\dot X\,Y+X\dot Y),
$$

$$
\ddot A
=\ddot A_0XY
+2\dot A_0(\dot X\,Y+X\dot Y)
+A_0(\ddot X\,Y+2\dot X\dot Y+X\ddot Y).
$$

这一个 product rule 分别作用于

| amplitude | $X$ | $Y$ |
|---|---|---|
| $A_{nn}=A_0N^2$ | $N$ | $N$ |
| $A_{\bar mn}=A_0\bar MN$ | $\bar M$ | $N$ |
| $A_{\bar m\bar m}=A_0\bar M^2$ | $\bar M$ | $\bar M$ |
| $\widetilde A_{nn}=A_0\widehat N^2$ | $\widehat N$ | $\widehat N$ |
| $\widetilde A_{\bar mn}=A_0\bar M\widehat N$ | $\bar M$ | $\widehat N$ |

$A_{ab}$ 是 A4 exact worldline amplitudes；带 tilde 的两个 amplitudes 是 compactified-$\Delta$ extension 所需的 horizon-regular amplitudes。

为使最终四块不依赖抽象 product-rule placeholders，实际使用的三个 extension amplitudes 及其 jets 显式为

$$
\begin{aligned}
\widetilde A_{nn}^{(0)}={}&A_0\widehat N^2,\\
\dot{\widetilde A}_{nn}={}&
\dot A_0\widehat N^2+2A_0\widehat N\dot{\widehat N},\\
\ddot{\widetilde A}_{nn}={}&
\ddot A_0\widehat N^2
+4\dot A_0\widehat N\dot{\widehat N}
+2A_0\left[(\dot{\widehat N})^2
+\widehat N\ddot{\widehat N}\right],
\end{aligned}
$$

$$
\begin{aligned}
\widetilde A_{\bar mn}^{(0)}={}&A_0\bar M\widehat N,\\
\dot{\widetilde A}_{\bar mn}={}&
\dot A_0\bar M\widehat N
+A_0\left(\dot{\bar M}\widehat N
+\bar M\dot{\widehat N}\right),\\
\ddot{\widetilde A}_{\bar mn}={}&
\ddot A_0\bar M\widehat N
+2\dot A_0\left(\dot{\bar M}\widehat N
+\bar M\dot{\widehat N}\right)\\
&+A_0\left(
\ddot{\bar M}\widehat N
+2\dot{\bar M}\dot{\widehat N}
+\bar M\ddot{\widehat N}\right),
\end{aligned}
$$

$$
\begin{aligned}
A_{\bar m\bar m}^{(0)}={}&A_0\bar M^2,\\
\dot A_{\bar m\bar m}={}&
\dot A_0\bar M^2+2A_0\bar M\dot{\bar M},\\
\ddot A_{\bar m\bar m}={}&
\ddot A_0\bar M^2
+4\dot A_0\bar M\dot{\bar M}
+2A_0\left[(\dot{\bar M})^2
+\bar M\ddot{\bar M}\right].
\end{aligned}
$$

右端的 $A_0$、$\widehat N,\bar M$ 及其 dot quantities 均已在本节前文显式定义。

### E.5.3 Azimuthal phase jets

按 A4 的 $e^{im\Phi}$ convention，定义

$$
\boxed{Q_{ab}=A_{ab}e^{-im\Phi_p}}.
$$

对 exact 或 tilde amplitude 使用同一公式：

$$
\boxed{
\dot Q_{ab}
=e^{-im\Phi_p}
(\dot A_{ab}-im\dot\Phi_pA_{ab})
},
$$

$$
\boxed{
\ddot Q_{ab}
=e^{-im\Phi_p}
\left[
\ddot A_{ab}
-2im\dot\Phi_p\dot A_{ab}
-im\ddot\Phi_pA_{ab}
-m^2\dot\Phi_p^2A_{ab}
\right]
}.
$$

后文记

$$
\widetilde Q_{nn}=\widetilde A_{nn}e^{-im\Phi_p}
=\frac{Q_{nn}}{h_p^2},
$$

$$
\widetilde Q_{\bar mn}
=\widetilde A_{\bar mn}e^{-im\Phi_p}
=\frac{Q_{\bar mn}}{h_p},
$$

$$
\widetilde Q_{\bar m\bar m}=Q_{\bar m\bar m}.
$$

在 $h_p=0$ 时，右侧 quotient 只表示 exterior identity；实际定义采用左侧 regular contractions，不能计算 $0/0$。

<a id="a5-extension"></a>

## E.6 Compactified-$\Delta$ off-worldline extension

### E.6.1 Field-point 与 worldline factors

选择以下 radial extension：

$$
\boxed{
F_{nn}^{\rm ext}
=\widehat\Delta(R)^2\widetilde Q_{nn}(T)G_RG_y
},
$$

$$
\boxed{
F_{\bar mn}^{\rm ext}
=\widehat\Delta(R)\widetilde Q_{\bar mn}(T)G_RG_y
},
$$

$$
\boxed{
F_{\bar m\bar m}^{\rm ext}
=Q_{\bar m\bar m}(T)G_RG_y
}.
$$

其中：

- field-point factors：$\rho(R,y),\bar\rho(R,y),\Delta(R),\Sigma(R,y)$、$j_T(R),j_R(R)$、$\widehat\Delta(R)$、$\widehat\Sigma(R,y)$、$c(y),\ell_s(y)$；
- worldline factors：$r_p,R_p,\Phi_p,u_p^A$ 及其 jets、$A_0,N_p,\widehat N_p,\bar M_p$、$Q_{ab},\widetilde Q_{ab}$；
- mixed kernel：$G_R(R-R_p)G_y(y-y_p)$，其 center 是 worldline data，argument 是 field point。

不能把 $\widehat\Delta(R)$ 偷换成 $\widehat\Delta_p$：前者规定 Gaussian 离开 worldline 后的 radial tetrad scaling，后者只属于 amplitude normalization。

### E.6.2 恢复 exact coefficient

在 $R=R_p$ 且 $h_p\ne0$，

$$
\widehat\Delta(R_p)^2\widetilde Q_{nn}=Q_{nn},
$$

$$
\widehat\Delta(R_p)\widetilde Q_{\bar mn}=Q_{\bar mn}.
$$

第三个 projection 显然不变。更强地，在 distribution level，

$$
\left[\frac{\widehat\Delta(R)}{h_p}\right]^k
Q(T)\delta(R-R_p)=Q(T)\delta(R-R_p),
\qquad k=1,2.
$$

因此在 Gaussian replacement 之前，这一 extension 与 A4 是同一个 spacetime distribution；对其再作任何有序 source derivatives 仍给出同一 distribution。horizon 点的等式由 $\widehat N$ 的连续定义取得，不通过 quotient 计算。

有限 $\sigma_R$ 时，不同但 worldline-equivalent 的 smooth extension 一般给出不同 source；这是真实的 regularization ambiguity，必须在 point-source limit 或 width variation 中量化，不能把历史代码选择当作唯一物理结论。

<a id="a5-gaussian-jets"></a>

## E.7 Moving Gaussian 的时间 jets

本节先保留一般 $y_p(T)$，随后取赤道 $y_p=0$。令

$$
x_R=R-R_p(T),\qquad x_y=y-y_p(T),
\qquad G=G_RG_y.
$$

固定 $\sigma_R,\sigma_y$ 时，

$$
\Gamma_1
=\frac{x_R\dot R_p}{\sigma_R^2}
+\frac{x_y\dot y_p}{\sigma_y^2},
$$

$$
\Gamma_2
=\Gamma_1^2
+\frac{x_R\ddot R_p-\dot R_p^2}{\sigma_R^2}
+\frac{x_y\ddot y_p-\dot y_p^2}{\sigma_y^2}.
$$

直接微分指数得到

$$
\dot G=\Gamma_1G,
\qquad
\ddot G=\Gamma_2G.
$$

对任意 phase-dressed amplitude $Z(T)$，定义

$$
f_Z^{(0)}=ZG,
$$

$$
\boxed{
f_Z^{(1)}=(\dot Z+Z\Gamma_1)G
},
$$

$$
\boxed{
f_Z^{(2)}=(\ddot Z+2\dot Z\Gamma_1+Z\Gamma_2)G
}.
$$

赤道 geodesic 有

$$
y_p=0,\qquad\dot y_p=0,\qquad\ddot y_p=0,
$$

故

$$
\Gamma_1=\frac{(R-R_p)\dot R_p}{\sigma_R^2},
$$

$$
\Gamma_2
=\frac{(R-R_p)^2\dot R_p^2}{\sigma_R^4}
+\frac{(R-R_p)\ddot R_p-\dot R_p^2}{\sigma_R^2}.
$$

若以后令 width 随时间或 grid 改变，上式不再成立；$\dot\sigma$、$\ddot\sigma$ 项必须显式加入。

对三个 projections，使用

$$
\widetilde f_{nn}^{(k)}
=f_{\widetilde Q_{nn}}^{(k)},
\qquad
\widetilde f_{\bar mn}^{(k)}
=f_{\widetilde Q_{\bar mn}}^{(k)},
$$

$$
f_{\bar m\bar m}^{(k)}
=f_{Q_{\bar m\bar m}}^{(k)},
\qquad k=0,1,2.
$$

由于 $\widehat\Delta(R)$ 是 time-independent field-point factor，完整 extended projection jets 是

$$
\partial_T^kF_{nn}^{\rm ext}
=\widehat\Delta(R)^2\widetilde f_{nn}^{(k)},
$$

$$
\partial_T^kF_{\bar mn}^{\rm ext}
=\widehat\Delta(R)\widetilde f_{\bar mn}^{(k)},
$$

$$
\partial_T^kF_{\bar m\bar m}^{\rm ext}
=f_{\bar m\bar m}^{(k)}.
$$

<a id="a5-four-blocks"></a>

## E.8 逐块消去所有 coordinate $\partial_T$

为避免把长 product rules 隐藏在未定义 placeholders 中，先定义只含 spatial derivatives 的 jet operators。对一列已知函数 $X^{(0)},X^{(1)},X^{(2)}$，令

$$
\mathfrak L_s[X]^{(k)}
=c\,\partial_yX^{(k)}
-iac\,X^{(k+1)}+\ell_sX^{(k)},
\qquad k=0,1,
$$

$$
\mathfrak J[X]^{(k)}
=j_TX^{(k+1)}+j_R\partial_RX^{(k)},
\qquad k=0,1.
$$

这不是新的物理算符，而是恒等式

$$
\widehat{\mathcal L}^{[m]}_sX^{(k)}
=\mathfrak L_s[X]^{(k)},
\qquad
\widehat{\mathcal J}_+X^{(k)}
=\mathfrak J[X]^{(k)}
$$

在 $\partial_TX^{(k)}=X^{(k+1)}$ 已解析代入后的记账形式。以下四块中只剩 $\partial_R,\partial_y$。

共同因子的来源可逐块追踪。使用

$$
\Delta^2=\frac{\widehat\Delta^2}{R^4},
\qquad
-\frac{16\pi\Sigma}{\Delta^2R}
=-\frac{16\pi R\widehat\Sigma}{\widehat\Delta^2},
$$

block 1 的两个 angular/time operators 不作用于 $\widehat\Delta(R)^2$；block 2 的 inner $\Delta^{-1}$ 与 projection 的一个 $\widehat\Delta(R)$ 相消；block 3 直接用其 outer $\Delta^2$；block 4 的 inner angular/time operator 先让 $\widehat\Delta(R)$ 穿过，再与 $\Delta^{-1}$ 相消。故每块的 $\widehat\Delta^2$ 都与 external prefactor 逐块相消，而不是在总和中相消。

### E.8.1 $nn$ / angular-angular block

定义

$$
U_1^{(k)}=\rho^{-2}\bar\rho^{-1}
\widetilde f_{nn}^{(k)},
\qquad k=0,1,2,
$$

$$
V_1^{(k)}=\mathfrak L_0[U_1]^{(k)},
\qquad k=0,1,
$$

$$
W_1^{(k)}=\rho^{-4}V_1^{(k)},
\qquad k=0,1.
$$

消去 $\widehat\Delta^2$ 后的 block kernel 是

$$
\boxed{
\mathfrak H_1
=-\rho^8\bar\rho\,
\mathfrak L_{-1}[W_1]^{(0)}
}.
$$

### E.8.2 $\bar mn$ / outer angular, inner radial block

原 inner factor

$$
\Delta^{-1}F_{\bar mn}^{\rm ext}
=\frac{R^2}{\widehat\Delta}
\widehat\Delta\widetilde Q_{\bar mn}G
=R^2\widetilde Q_{\bar mn}G
$$

先解析相消。定义

$$
U_2^{(k)}=\rho^{-2}\bar\rho^{-2}R^2
\widetilde f_{\bar mn}^{(k)},
$$

$$
V_2^{(k)}=\mathfrak J[U_2]^{(k)},
\qquad
W_2^{(k)}=\rho^{-4}\bar\rho^2V_2^{(k)},
$$

其中 $k=0,1$；$U_2$ 另需 $k=2$ 以构造 $V_2^{(1)}$。则

$$
\boxed{
\mathfrak H_2
=-\frac{\rho^8\bar\rho}{\sqrt2R^4}
\mathfrak L_{-1}[W_2]^{(0)}
}.
$$

### E.8.3 $\bar m\bar m$ / radial-radial block

定义

$$
U_3^{(k)}=\rho^{-2}\bar\rho
f_{\bar m\bar m}^{(k)},
\qquad k=0,1,2,
$$

$$
V_3^{(k)}=\mathfrak J[U_3]^{(k)},
\qquad k=0,1,
$$

$$
W_3^{(k)}=\rho^{-4}V_3^{(k)},
\qquad k=0,1.
$$

于是

$$
\boxed{
\mathfrak H_3
=-\frac{\rho^8\bar\rho}{2R^4}
\mathfrak J[W_3]^{(0)}
}.
$$

### E.8.4 $\bar mn$ / outer radial, inner angular block

内层 angular operator 不作用于 field-point-only 的 $\widehat\Delta(R)$，故该因子穿过内层后与 $\Delta^{-1}$ 相消。定义

$$
U_4^{(k)}=\rho^{-2}\bar\rho^{-2}
\widetilde f_{\bar mn}^{(k)},
\qquad k=0,1,2,
$$

$$
V_4^{(k)}=\mathfrak L_{-1}[U_4]^{(k)},
\qquad k=0,1,
$$

$$
W_4^{(k)}=\rho^{-4}\bar\rho^2R^2V_4^{(k)},
\qquad k=0,1.
$$

于是

$$
\boxed{
\mathfrak H_4
=-\frac{\rho^8\bar\rho}{\sqrt2R^4}
\mathfrak J[W_4]^{(0)}
}.
$$

### E.8.5 四块完整数值源

四块共享已经解析正则化的 field-point prefactor：

$$
\boxed{
S_{m,A}^{\rm num}
\equiv S_{m,A}^{G,\rm ext}
=-16\pi R\widehat\Sigma(R,y)\,\mathfrak H_A,
\qquad A=1,2,3,4
}.
$$

因此

$$
\boxed{
S_m^{\rm num}
=\sum_{A=1}^4S_{m,A}^{\rm num}
}.
$$

这些表达式不含 Dirac delta，不含 coordinate $\partial_T$，也不含未定义 amplitude jets；仅保留按原 nesting 作用于已知 Gaussian products 的 $\partial_R$、$\partial_y$。它们在 $0<R<R_H$、$|y|<1$ 上是 smooth functions。这里还没有把 spatial derivatives 离散化或展开成 stencil。

## E.9 每块独立的 horizon / SCRI+ radial regularity

### E.9.1 Future horizon

假设非极端 Kerr，$R_H=L^2/r_+>0$。在 $R\to R_H^-$：

- $R^{-4}$ 有限；
- $\rho,\bar\rho,j_T,j_R,\widehat\Sigma$ 有限；
- $\widetilde Q_{nn},\widetilde Q_{\bar mn}$ 由 $\widehat N$ 直接定义，有限；
- $G$ 及其 time/spatial derivatives 有限；
- 所有 $\widehat\Delta^{-1}$ 已在 E.8 逐块解析消去。

故每个 $\mathfrak H_A$ 和 $S_{m,A}^{\rm num}$ 单独具有普通 Taylor/Laurent expansion，$\widehat\Delta$ 的负幂系数逐块为零。不需要也不允许用 $\sum_A S_A$ 的块间相消建立 regularity。

### E.9.2 SCRI+

在固定 $|y|<1$、$R\to0^+$，

$$
\rho=O(R),\quad\bar\rho=O(R),
\quad j_T=O(1),\quad j_R=O(R^2),
$$

$$
\widehat\Sigma=O(1),
\qquad -16\pi R\widehat\Sigma=O(R).
$$

若 Gaussian 和 worldline amplitudes 在该 field-point limit 中取 $O(1)$，则逐块最坏幂次为

| block | $U_A$ | $V_A$ | $W_A$ | $\mathfrak H_A$ | $S_{m,A}^{\rm num}$ |
|---|---:|---:|---:|---:|---:|
| 1: $nn$ | $O(R^{-3})$ | $O(R^{-3})$ | $O(R^{-7})$ | $O(R^2)$ | $O(R^3)$ |
| 2: $\bar mn$ LJ | $O(R^{-2})$ | $O(R^{-2})$ | $O(R^{-4})$ | $O(R)$ | $O(R^2)$ |
| 3: $\bar m\bar m$ JJ | $O(R^{-1})$ | $O(R^{-1})$ | $O(R^{-5})$ | $O(1)$ | $O(R)$ |
| 4: $\bar mn$ JL | $O(R^{-4})$ | $O(R^{-4})$ | $O(R^{-4})$ | $O(R)$ | $O(R^2)$ |

$j_R\partial_R$ 最多把 operand 的幂降低一阶后再乘 $R^2$，因此不比表中由 $j_T\partial_T$ 给出的 leading order 更坏。表中每一行的 final source 都没有 $R$ 的负幂，并分别趋于零；这是一项逐块 analytic power-counting proof，不是对 total source 的经验 cancellation。

该结论只处理 radial boundaries。$c^{-1}$ 在 $y=\pm1$ 的问题不包含在此 Laurent 证明中。

<a id="a5-contract"></a>

## E.10 A5 的精确输入、输出与逐块判据

A5 的 frozen theory input 是 [附录 D.5](W02_1_4_FIXED_M_Y_DERIVATION.md#a4-y-source) 的四个 exact blocks，加上本附录明确选择的 fixed-width Gaussian 与 compactified-$\Delta$ extension。每个 stage 的 independent scalar inputs 是

$$
\boxed{
\{M,a,L,\mu,m,E,L_z,r_p,\Phi_p,\sigma_R,\sigma_y\}
},
$$

并取 field-point variables $(R,y)$。$R_p,u_p^A$、全部 radial/time jets、$\widehat N_p,\bar M_p$、amplitude/phase/Gaussian jets 都必须由这些 scalars 和本文 formulas 生成，不能作为额外未验证输入。

A5 的逐块 theory outputs 是

$$
\boxed{
\left\{
S_{m,1}^{\rm num},
S_{m,2}^{\rm num},
S_{m,3}^{\rm num},
S_{m,4}^{\rm num}
\right\}
},
$$

连同中间 projection-jet triples

$$
\left\{
\widetilde f_{nn}^{(0,1,2)},
\widetilde f_{\bar mn}^{(0,1,2)},
f_{\bar m\bar m}^{(0,1,2)}
\right\}.
$$

对每个 $A=1,\ldots,4$，验证路线必须分别：

1. 从对应 A4 exact block 代入 extension 后，让 Mathematica 直接执行全部 coordinate $\partial_T$；
2. 独立用本文 worldline/contraction/amplitude/phase/Gaussian jets 构造 $\mathfrak H_A$；
3. 比较完整 $S_{m,A}^{\rm num}$，要求 residual 严格为 $0$。

每块 output 必须满足：不含 Dirac delta、不含未求值 coordinate $\partial_T$、不含未定义 dot quantities；只允许 $\partial_R,\partial_y$ 按 A4 ordering 作用于已知 smooth operands。四块通过后才比较

$$
S_m^{\rm num}=\sum_{A=1}^4S_{m,A}^{\rm num}.
$$

Finite-domain normalization、axis endpoint completion、width convergence 与 turning-point alternative parameterization 不属于上述 time-algebra residual，必须作为独立边界保留。

实际 A5 staged validation 从上述 inputs 出发完成了 common jets 与四个 blocks 的独立比较：所有 common-jet、worldline-recovery、$S_{m,A}^{\rm num}$ 和 total-source residual 均为 $0$，且生成路线未调用理论 $\mathfrak L_s,\mathfrak J,\mathfrak H_A$ helpers。机器汇总见 [A5 result](../../../../data/kerr_point_particle_verification/results/a5/a5_endpoint_summary.json)，给出 `A5FullEndpointVerifiedWithinInteriorScope=True`。

## E.11 历史实现的使用与证据强度

为识别 repository plan 所称的“当前 compactified-$\Delta$ extension”，本文在 C01 已授权范围内只读取了目标仓库

- `NewPointParticlesm2Evolution/appendix/compactified_delta_gaussian_extension.md`；
- `NewPointParticlesm2Evolution/appendix/fixed_m_gaussian_source_after_t_derivatives.md`；
- `src/sminus2_point_particle/source.py`、`trajectory.py`、`geometry.py`；
- 与上述 source 直接对应的 generated Wolfram Language 和历史 boundary/result records。

历史材料提示了 $\widehat\Delta^2,\widehat\Delta,1$ 三种 projection extension、regular $\widehat N$ 和 time-jet nesting。本文从 A4 exact source 重新推得 E.6--E.9 的 distribution identity、逐块 cancellation 与 power counting。历史 JSON 所记的 zero residual 和 deterministic Laurent check 没有在本仓库重放，本附录不声称执行或验证了它们；尤其历史 boundary check 对 total source 的一个确定 stage 取值，不能替代本仓库已经执行的 A5 common-jets 与四块 interface checks，也不能替代 W02.3 的 production boundary/convergence tests。

## E.12 开放问题与可支持结论

在 A4、A1--A3 条件成立且 $\mathcal V>0$、$P+\mathcal V\ne0$、$u^T>0$、固定 widths、$0<R<R_H$、$|y|<1$ 的范围内，本文支持：

1. compactified-$\Delta$ radial extension 在 distribution level 与 A4 exact source 等价；
2. 它在 $R=R_p$ 恢复三个 exact projection coefficients；
3. $R_p,\Phi_p$、contractions、amplitudes、phase 和 moving Gaussian 的 0/1/2 阶 time jets 已显式闭合；
4. 四个 nested source blocks 的 coordinate $\partial_T$ 已全部消去，并定义为 $S_{m,A}^{\rm num}$；
5. 四块分别在 horizon 与 SCRI+ 无 radial Laurent 负幂。

以下问题仍开放，可能改变有限宽 source 或物理结果：

- **A4 gate：** `A4FullEndpointVerified=True`；A5 继承该 exact source、A1 边界和 local assumptions。
- **A1 convention gate：** metric/Weyl-scalar bridge、sourced units 与 Ripley operator 独立转录尚未全部关闭；整体 source sign/normalization 仍继承该边界。
- **Finite-domain normalization：** full-line Gaussian 在有限 radial/angular domain 上不是严格单位积分；domain renormalization 会增加 time derivatives，不能在数值实现中暗中加入。
- **Angular-axis extension：** 对 $m=2,4$，普通 $G_y$ 的 axis tail 不会自行抵消 $\ell_s\sim c^{-1}$；而 A4 的 field axis powers 并不足以唯一规定三个 source projections 的 spin-frame completion。因此目前不能宣称 $S_m^{\rm num}$ 在 $y=\pm1$ 光滑。必须另外推导一个与 exact delta 等价、与各 projection spin weight 相容的 angular off-worldline extension，或明确采用不含 axis endpoints 的 angular representation；该选择属于数值接口的独立 closure。
- **Width ambiguity：** compactified-$\Delta$ 只是一个 worldline-equivalent extension；有限 $\sigma_R$ 的结果必须通过 width dependence / point-source limit 评估。
- **Turning point / branch：** 本文 radial derivative formulas 在 $\mathcal V=0$ 处奇异；若 W02.3 的初值落在精确 turning point，必须换用适合该点的参数化，不能直接使用本附录 jets。

因此，本附录给出的 interior $|y|<1$ A5 time-algebra target 已通过逐块验证；但在 finite-domain normalization、angular-axis interface、turning-point parameterization 与 width convergence 分别关闭前，仍不应作为完整 theory-to-code contract。
