# 附录 C：赤道 Kerr 点粒子、constant-$T$ stress tensor 与四块 source

- 主文：[DERIVATION.md](../DERIVATION.md)
- 范围：W02.1.3
- 任务：T006
- 状态：A3 local algebraic endpoint 已验证；继承 A1；global worldline/turning-point 边界开放
- Primary sources：Sasaki--Tagoshi `P08` Eqs. (26)--(32)；Mino--Brink `P09` Eqs. (2.1)--(2.3)
- Inherited input：A2 的完整 source functional 与 sourced Ripley equation，见 [附录 B](SASAKI_SOURCE.md)

## C.1 证据层级与 A3 起点

本附录区分三类内容：

- `primary`：`P08` 的 point-particle stress tensor、Kerr geodesic first integrals 与 $T_{nn},T_{\bar mn},T_{\bar m\bar m}$，以及 `P09` 的 Kerr plunge geodesic equations；
- `derived`：赤道化简、hyperboloidal coordinate transformation、constant-$T$ delta integration、determinant/delta Jacobians、horizon-regular rewriting 和具体 point-particle four blocks；
- `inherited`：A2 已建立的四个外部 source operators 以及条件 sourced equation

$$
\widehat T_H
\left[T_{nn}^H,T_{\bar mn}^H,T_{\bar m\bar m}^H\right]
=\sum_{A=1}^4\widehat T_{H,A},
$$

$$
\mathcal O_{\rm R}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}\widehat T_H.
$$

A3 不重新推导 A2 operators。其任务是：给定 W02.3 将提供的 $(E,L_z,r_0)$，从赤道 timelike geodesic 和 invariant point-particle stress tensor 独立得到 $T_{nn}^{H,\rm pp},T_{\bar mn}^{H,\rm pp},T_{\bar m\bar m}^{H,\rm pp}$，再逐项代入四块 A2 functional。

采用 $(-,+,+,+)$ Sasaki convention，$E$ 和 $L_z$ 是除以粒子质量 $\mu$ 后的 specific constants；Carter constant 取 $Q=0$。背景量为

$$
\Sigma=r^2+a^2\cos^2\theta,
\qquad
\Delta=r^2-2Mr+a^2.
$$

## C.2 赤道 Kerr timelike geodesic

### C.2.1 BL first integrals

定义

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

`P08` Eqs. (27)--(28) 与 `P09` Eqs. (2.1)--(2.3) 给出 Kerr geodesic first integrals。取

$$
\theta_p=\frac\pi2,
\qquad Q=0,
\qquad \Sigma_p=r_p^2,
$$

以及 future-directed inward branch 后，完整 BL four-velocity 为

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

这些 expressions 满足 timelike normalization 与 conserved quantities

$$
g_{\mu\nu}u^\mu u^\nu=-1,
\qquad E=-u_t,
\qquad L_z=u_\phi.
$$

适用域要求 $\mathcal R\ge0$。负号选定 inward branch；在精确 turning point $\mathcal R=0$ 上 velocity 连续，但以后若需要 $r$ derivatives，必须另行处理 $\sqrt{\mathcal R}$ 的可微性。外部 BL expressions 取 $r>r_+$、$\Delta>0$；horizon limit 用下一节的 regular form 定义。

Future direction 要求 $u^t>0$。对穿过 future horizon 的 branch，另要求

$$
P_+\equiv P(r_+)
=(r_+^2+a^2)(E-\Omega_HL_z)>0,
\qquad
\Omega_H=\frac{a}{r_+^2+a^2}.
$$

### C.2.2 Hyperboloidal velocity 与 horizon-regular form

坐标变换为

$$
T=t+h(r),\qquad R=\frac{L^2}{r},
\qquad \Phi=\phi+q(r),\qquad y=-\cos\theta,
$$

$$
h'(r)=\frac{r^2+a^2}{\Delta}-2-\frac{4M}{r},
\qquad q'(r)=\frac a\Delta.
$$

Vector transformation 给出

$$
\boxed{
u^T=u^t+h'u^r,
\quad u^R=-\frac{L^2}{r^2}u^r,
\quad u^\Phi=u^\phi+\frac a\Delta u^r,
\quad u^y=\sin\theta\,u^\theta=0
}.
$$

由

$$
(P-\sqrt{\mathcal R})(P+\sqrt{\mathcal R})
=P^2-\mathcal R=\Delta C
$$

定义

$$
\boxed{
K(r)=\frac{P-\sqrt{\mathcal R}}{\Delta}
=\frac{C}{P+\sqrt{\mathcal R}}
}.
$$

第二式在 $P+\sqrt{\mathcal R}\ne0$ 时是 regular definition。于是

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
\qquad u^\Phi=\frac{B+aK}{r^2},
\qquad u^y=0
}.
$$

在 future-horizon branch 上 $\sqrt{\mathcal R}\to P_+>0$，故

$$
K\longrightarrow\frac{C(r_+)}{2P_+}
$$

且 $(u^T,u^R,u^\Phi)$ 均有限。W02.3 应直接使用这一 regular form，不能分别计算 divergent BL terms 后依赖 floating-point cancellation。

给定 $(E,L_z,r_0)$ 后，$T$-parameterized worldline 由

$$
\frac{dr_p}{dT}=\frac{u_p^r}{u_p^T},
\qquad
\frac{dR_p}{dT}=\frac{u_p^R}{u_p^T},
\qquad
\frac{d\Phi_p}{dT}=\frac{u_p^\Phi}{u_p^T}
$$

和 $r_p(T_0)=r_0$ 唯一指定；$R_p(T_0)=L^2/r_0$，$\Phi_p(T_0)$ 只固定轴对称背景中的 azimuthal origin。这里没有选择 $E,L_z,r_0$ 的 ISCO/transition prescription，该输入仍属于 W02.3。

## C.3 Invariant point-particle stress tensor

### C.3.1 Constant-$T$ integration 与 absolute value

Invariant starting point 是

$$
\begin{aligned}
T^{AB}(X)={}&\mu\int
\frac{u^Au^B}{\sqrt{-g_H}}
\delta\!\left(T-T_p(\tau)\right)
\delta\!\left(R-R_p(\tau)\right)\\
&\times\delta\!\left(\theta-\theta_p(\tau)\right)
\delta\!\left(\Phi-\Phi_p(\tau)\right)d\tau.
\end{aligned}
$$

若 constant-$T$ hypersurface 与 worldline 只有一个 transverse intersection $\tau_p(T)$，则

$$
\int d\tau\,\delta\!\left(T-T_p(\tau)\right)
=\frac1{|u_p^T|}.
$$

所以

$$
\boxed{
T_H^{AB}
=\frac{\mu u_p^Au_p^B}{\sqrt{-g_H}\,|u_p^T|}
\delta(R-R_p)\delta(\theta-\theta_p)\delta(\Phi-\Phi_p)
}.
$$

本文采用 future-directed branch $u_p^T>0$，于是 $|u_p^T|=u_p^T$。若 $u^T=0$，slice 不再 transverse，以上 reduction 失效。

### C.3.2 Determinant、coordinate measure 与 BL delta route

BL determinant 与 inverse-coordinate Jacobian 分别为

$$
\det g_{tr\theta\phi}=-\Sigma^2\sin^2\theta,
$$

$$
\left|\det\frac{\partial(t,r,\theta,\phi)}
{\partial(T,R,\theta,\Phi)}\right|
=\left|\frac{dr}{dR}\right|=\frac{r^2}{L^2}.
$$

因此

$$
\boxed{
\sqrt{-g_H}=\frac{\Sigma r^2\sin\theta}{L^2}
}.
$$

三个 deltas 相对于 coordinate measure $dR\,d\theta\,d\Phi$ 定义。

作为独立 BL-to-hyperboloidal route，在 fixed $T$ 上令

$$
F(r)=r-r_p\!\left(T-h(r)\right).
$$

在 root $r=r_p$，

$$
F'(r_p)=1+h'_p\frac{u_p^r}{u_p^t}
=\frac{u_p^T}{u_p^t}.
$$

对 $u_p^t>0,u_p^T>0$，再用 $|dR/dr|=L^2/r_p^2$，得到

$$
\boxed{
\frac1{u_p^t}\delta\!\left(r-r_p(t)\right)
=\frac{L^2}{r_p^2u_p^T}
\delta\!\left(R-R_p(T)\right)
}.
$$

在 radial support 上，

$$
\Phi-q(r)-\phi_p\!\left(T-h(r)\right)
=\Phi-\Phi_p(T),
$$

且 argument 对 $\Phi$ 的 derivative 为 $1$，故 azimuthal delta 无额外 Jacobian；$\theta$ 在本阶段未变换。这条 route 与 invariant constant-$T$ integration 给出同一 distribution weight。

## C.4 Kinnersley covectors 与 contractions

### C.4.1 Covector pullback

从 `P08` Eqs. (4)--(5) 的 Kinnersley vectors 在 BL metric 中 lowering，再按 covector rule pull back。分量顺序为 $(T,R,\theta,\Phi)$：

$$
\boxed{
n_A^{\rm K}
=\left(
-\frac\Delta{2\Sigma},
\frac{\Delta r^2}{L^2\Sigma}\left(1+\frac{2M}{r}\right),
0,
\frac{a\Delta\sin^2\theta}{2\Sigma}
\right)
},
$$

$$
\boxed{
\bar m_A^{\rm K}
=\frac1{\sqrt2(r-ia\cos\theta)}
\left(
ia\sin\theta,
-\frac{iar^2\sin\theta}{L^2}\left(2+\frac{4M}{r}\right),
\Sigma,
-i(r^2+a^2)\sin\theta
\right)
}.
$$

这里只变换同一 Kinnersley tetrad 的 coordinate components；没有额外应用 Ripley regular-tetrad boost/spin rotation。定义 scalars

$$
N_p=n_A^{\rm K}u_p^A,
\qquad \bar M_p=\bar m_A^{\rm K}u_p^A.
$$

显式地，

$$
\begin{aligned}
N_p={}&-\frac{\Delta_p}{2\Sigma_p}u_p^T
+\frac{\Delta_pr_p^2}{L^2\Sigma_p}
\left(1+\frac{2M}{r_p}\right)u_p^R\\
&+\frac{a\Delta_p\sin^2\theta_p}{2\Sigma_p}u_p^\Phi,
\end{aligned}
$$

$$
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
$$

### C.4.2 赤道 $N_p,\bar M_p$

代入 geodesic velocities 得到

$$
\boxed{
N_p=-\frac{P_p-\sqrt{\mathcal R_p}}{2r_p^2}
=-\frac{\Delta_pK_p}{2r_p^2}
},
$$

$$
\boxed{
\bar M_p=\frac{i(aE-L_z)}{\sqrt2\,r_p}
=-\frac{iB}{\sqrt2\,r_p}
}.
$$

$N_p=O(\Delta_p)$ 对 future-horizon regularity 至关重要。后续 radial extension 可定义

$$
\widehat N_p\equiv\frac{N_p}{R_p^2\Delta_p}
=-\frac{K_p}{2L^4},
$$

但 exact A3 projections 仍用 $N_p$。

## C.5 $T_{nn}^{H,\rm pp},T_{\bar mn}^{H,\rm pp},T_{\bar m\bar m}^{H,\rm pp}$

定义

$$
\delta_p^{(3)}
\equiv\delta(R-R_p(T))\delta(\theta-\theta_p)
\delta(\Phi-\Phi_p(T)),
$$

$$
\boxed{
\mathcal W_p^{(\theta)}
=\frac{\mu L^2}
{\Sigma_pr_p^2\sin\theta_p\,u_p^T}\delta_p^{(3)}
}.
$$

三个 Kinnersley 分量是

$$
\boxed{
T_{nn}^{H,\rm pp}=\mathcal W_p^{(\theta)}N_p^2,
\qquad
T_{\bar mn}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_pN_p,
\qquad
T_{\bar m\bar m}^{H,\rm pp}=\mathcal W_p^{(\theta)}\bar M_p^2
}.
$$

对 equatorial orbit，

$$
\mathcal W_p^{(\theta)}
=\frac{\mu L^2}{r_p^4u_p^T}
\delta(R-R_p)\delta\!\left(\theta-\frac\pi2\right)
\delta(\Phi-\Phi_p).
$$

## C.6 逐项代入后的 point-particle four blocks

以下把 projections 逐项代入 A2 的四个有序 blocks。为使 endpoint 可独立录入，保留最小 operator definitions：

$$
\rho=(r-ia\cos\theta)^{-1},
\qquad
\bar\rho=(r+ia\cos\theta)^{-1},
\qquad r=\frac{L^2}{R},
$$

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

所有 $\rho,\bar\rho,\Delta$ 是 field-point functions；$\mathcal W_p^{(\theta)},N_p,\bar M_p$ 含完整 worldline amplitudes 与 distribution support。Operators 由内向外作用于其右侧全部 smooth factors 和 distributions。

$$
\boxed{
\widehat T_{H,1}^{\rm pp}
=-\rho^8\bar\rho\,
\widehat{\mathcal L}_{-1}
\left[
\rho^{-4}\widehat{\mathcal L}_0
\left(\rho^{-2}\bar\rho^{-1}
\mathcal W_p^{(\theta)}N_p^2\right)
\right]
}.
$$

$$
\boxed{
\begin{aligned}
\widehat T_{H,2}^{\rm pp}
={}&-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal L}_{-1}
\Bigg[
\rho^{-4}\bar\rho^2\widehat{\mathcal J}_+
\Bigg(\rho^{-2}\bar\rho^{-2}\Delta^{-1}
\mathcal W_p^{(\theta)}\bar M_pN_p\Bigg)
\Bigg].
\end{aligned}
}
$$

$$
\boxed{
\begin{aligned}
\widehat T_{H,3}^{\rm pp}
={}&-\frac12\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\Bigg[
\rho^{-4}\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho\,
\mathcal W_p^{(\theta)}\bar M_p^2\right)
\Bigg].
\end{aligned}
}
$$

$$
\boxed{
\begin{aligned}
\widehat T_{H,4}^{\rm pp}
={}&-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\Bigg[
\rho^{-4}\bar\rho^2\Delta^{-1}
\widehat{\mathcal L}_{-1}
\left(\rho^{-2}\bar\rho^{-2}
\mathcal W_p^{(\theta)}\bar M_pN_p\right)
\Bigg].
\end{aligned}
}
$$

两个 mixed blocks 的 ordering 不同，不能交换。定义

$$
\widehat T_H^{\rm pp}=\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}.
$$

完整 unprojected point-particle sourced equation 是

$$
\boxed{
\mathcal O_{\rm R}[\psi_4]=S_{\rm R}^{\rm pp},
\qquad
S_{\rm R}^{\rm pp}
=-\frac{16\pi\Sigma}{\Delta^2R}
\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}
}.
$$

外部 $\Sigma/(\Delta^2R)$ 也是 field-point multiplier。该 equation 仍继承 A1 的 signature/Weyl-scalar、sourced-unit 和 Ripley-PDF transcription 边界。

## C.7 A3 closure 结果与未决边界

A3 Mathematica check 已从两条独立路线计算同一 endpoint：

1. 从 invariant four-dimensional stress tensor 做 constant-$T$ root integration、determinant/Jacobian transformation 和 tetrad contractions；
2. 从 `P08` Eq. (26) 的 BL time-sliced stress tensor，对 fixed-$T$ moving radial root、$R=L^2/r$、azimuthal delta 和 covectors 作完整 pullback。

比较对象包括 BL/hyperboloidal four-velocity、$\sqrt{-g_H}$、全部 spatial deltas、$N_p,\bar M_p$、$T_{nn}^{H,\rm pp},T_{\bar mn}^{H,\rm pp},T_{\bar m\bar m}^{H,\rm pp}$，以及 $N_p,\bar M_p$ 的赤道公式。两条路线和 tetrad-projection stage 的 residual 均为 $0$；summary 给出 `A3FullEndpointVerified=True`。A3 在这三个 Kinnersley 分量处结束；四个 point-particle blocks 与 $S_{\rm R}^{\rm pp}$ 由它们和已验证的 A2 source functional 组合得到，不在 A3 重验。执行证据见 [A3 endpoint summary](../../../../data/kerr_point_particle_verification/results/a3/a3_endpoint_summary.json)。

已保存 residual 使用的 assumptions 为：

- $M>0,L>0,|a|<M$；
- exterior derivation 中 $r>r_+$、$\Delta>0$，horizon 用单侧 regular limit；
- $\mathcal R\ge0$，inward sign 已固定；
- $u^t>0,u^T>0$，constant-$T$ intersection unique and transverse；
- $P+\sqrt{\mathcal R}\ne0$；future-horizon crossing 取 $P_+>0$；
- equatorial timelike geodesic，$Q=0$，$u^\theta=u^y=0$。

该结论是上述 assumptions 下的 local algebraic closure，不证明 global worldline monotonicity、global unique slice intersection 或 turning-point differentiability，并继续继承 A1 convention/tetrad boundary。
