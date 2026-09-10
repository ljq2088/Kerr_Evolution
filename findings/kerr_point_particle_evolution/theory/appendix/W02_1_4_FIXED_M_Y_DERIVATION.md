# 附录 D：fixed-$m$ 与 $y=-\cos\theta$ exact distributional equation

- 主文：[DERIVATION.md](../DERIVATION.md)
- 范围：W02.1.4
- 任务：T006
- 状态：A4 完整 $\theta$ stage 与未投影点粒子方程到 $y$ 终点的 residual 均为 $0$；继承 A1 boundary
- 起点：[Section V frozen unprojected interface](../DERIVATION.md#v-完整未分解-point-particle-equation冻结接口)
- 上游：[A2 source](SASAKI_SOURCE.md)、[A3 point-particle components](POINT_PARTICLE_PROJECTIONS.md)
- 排除：Gaussian、off-worldline extension、time jets、spatial derivative expansion 和任何离散化

## D.1 起点、对象与工作边界

起点是已经冻结的完整 $(T,R,\theta,\Phi)$ 方程

$$
\boxed{
\mathcal O_{\rm R}^{(-2)}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}
\sum_{A=1}^4\widehat T_{H,A}^{\rm pp}
}.
$$

其中 $\mathcal O_{\rm R}^{(-2)}$ 是 [附录 A](RIPLEY_OPERATOR.md#ripley-eq-22) 的完整 Ripley operator 取 $s=-2$；四个 point-particle blocks 是 [附录 C.6](POINT_PARTICLE_PROJECTIONS.md#c6-逐项代入后的-point-particle-four-blocks) 的完整 expressions。

本附录只作两个连续变换：

$$
(T,R,\theta,\Phi)
\longrightarrow(T,R,\theta;m)
\longrightarrow(T,R,y;m),
\qquad y=-\cos\theta.
$$

四块 source 的 outer-to-inner ordering 始终不变。$\partial_R$、$\partial_\theta$ 或 $\partial_y$ 保持 nested operator form，不展开为大量 delta derivatives。

记

$$
r=\frac{L^2}{R},
\qquad
\Delta=r^2-2Mr+a^2.
$$

## D.2 第一步：fixed-$m$ Fourier projection

### D.2.1 Fourier pair 与 periodic delta

采用

$$
\boxed{
f(T,R,\theta,\Phi)
=\sum_{m\in\mathbb Z}f_m(T,R,\theta)e^{im\Phi}
},
$$

$$
\boxed{
f_m(T,R,\theta)
=\frac1{2\pi}\int_0^{2\pi}
f(T,R,\theta,\Phi)e^{-im\Phi}\,d\Phi
}.
$$

Kerr background 和所有 unprojected coefficients 与 $\Phi$ 无关，因此 mode projection 与 $\partial_T,\partial_R,\partial_\theta$ 交换，并有

$$
\boxed{
\partial_\Phi\longrightarrow im,
\qquad
\partial_\Phi^2\longrightarrow-m^2
}.
$$

Point-particle azimuthal delta 必须使用同一 Fourier convention：

$$
\boxed{
\delta_{2\pi}(\Phi-\Phi_p)
=\frac1{2\pi}\sum_{m\in\mathbb Z}
e^{im\Phi}e^{-im\Phi_p}
}.
$$

故每个 source mode 含且只含一次

$$
\boxed{
\frac{e^{-im\Phi_p(T)}}{2\pi}
}.
$$

### D.2.2 $(T,R,\theta;m)$ operators

定义 fixed-$m$ angular field operator

$$
\boxed{
\mathcal A_{s,m}^{(\theta)}[f]
=\frac1{\sin\theta}\partial_\theta
\left(\sin\theta\,\partial_\theta f\right)
+\left[
s-\frac{(m+s\cos\theta)^2}{\sin^2\theta}
\right]f
}.
$$

Source operators 变为

$$
\boxed{
\widehat{\mathcal L}_s^{[m,\theta]}
=\partial_\theta+m\csc\theta
-ia\sin\theta\,\partial_T+s\cot\theta
},
$$

$$
\boxed{
\widehat{\mathcal J}_+
=-\left(2+\frac{4M}{r}\right)\partial_T
-\frac{L^2}{r^2}\partial_R
}.
$$

$\widehat{\mathcal J}_+$ 的 $\partial_\Phi$ 已在 A2 coordinate pullback 中相消，因此 fixed-$m$ 后没有新增 $m$ term。

<a id="a4-theta-operator"></a>

### D.2.3 完整 $(T,R,\theta;m)$ field operator

对任意 spin weight $s$，完整 fixed-$m$ operator 是

$$
\begin{aligned}
\mathcal O_{s,m}^{(T,R,\theta)}[\psi_m]={}&
\left[
8M\left(2M-\frac{a^2R}{L^2}\right)
\left(1+\frac{2MR}{L^2}\right)-a^2\sin^2\theta
\right]\partial_T^2\psi_m
\\
&-2\left[
L^2-(8M^2-a^2)\frac{R^2}{L^2}
+\frac{4a^2MR^3}{L^4}
\right]\partial_T\partial_R\psi_m
\\
&-\left(L^2-2MR+\frac{a^2R^2}{L^2}\right)
\frac{R^2}{L^2}\partial_R^2\psi_m
-\mathcal A_{s,m}^{(\theta)}[\psi_m]
\\
&+2iam\left(1+\frac{4MR}{L^2}\right)\partial_T\psi_m
+\frac{2iamR^2}{L^2}\partial_R\psi_m
\\
&+2\left[
2M\left(-s+2(s+2)\frac{MR}{L^2}
-3\frac{a^2R^2}{L^4}\right)
-\frac{a^2R}{L^2}+isa\cos\theta
\right]\partial_T\psi_m
\\
&+2R\left[-(1+s)+(s+3)\frac{MR}{L^2}
-2\frac{a^2R^2}{L^4}\right]\partial_R\psi_m
\\
&+\frac{2iamR}{L^2}\psi_m
+2\left[(1+s)\frac{MR}{L^2}
-\frac{a^2R^2}{L^4}\right]\psi_m.
\end{aligned}
$$

W02 使用 $s=-2$。

### D.2.4 Fixed-$m$ point-particle components

在 $dR\,d\theta$ measure 下定义

$$
\boxed{
\mathcal W_m^{(\theta)}
=\frac{\mu L^2}
{2\pi\Sigma_pr_p^2\sin\theta_p\,u_p^T}
e^{-im\Phi_p(T)}
\delta(R-R_p(T))\delta(\theta-\theta_p)
}.
$$

于是

$$
\boxed{
T_{nn,m}^{H,\theta}=\mathcal W_m^{(\theta)}N_p^2,
\qquad
T_{\bar mn,m}^{H,\theta}=\mathcal W_m^{(\theta)}\bar M_pN_p,
\qquad
T_{\bar m\bar m,m}^{H,\theta}=\mathcal W_m^{(\theta)}\bar M_p^2
}.
$$

<a id="a4-theta-blocks"></a>

### D.2.5 四个 $(T,R,\theta;m)$ exact blocks

这里

$$
\rho=(r-ia\cos\theta)^{-1},
\qquad
\bar\rho=(r+ia\cos\theta)^{-1}.
$$

保持 A2 ordering：

$$
\begin{aligned}
\widehat T_{m,1}^{(\theta),\rm dist}={}&
-\rho^8\bar\rho\,
\widehat{\mathcal L}_{-1}^{[m,\theta]}
\left[
\rho^{-4}\widehat{\mathcal L}_0^{[m,\theta]}
\left(\rho^{-2}\bar\rho^{-1}T_{nn,m}^{H,\theta}\right)
\right],
\\[1ex]
\widehat T_{m,2}^{(\theta),\rm dist}={}&
-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal L}_{-1}^{[m,\theta]}
\left[
\rho^{-4}\bar\rho^2\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho^{-2}\Delta^{-1}
T_{\bar mn,m}^{H,\theta}\right)
\right],
\\[1ex]
\widehat T_{m,3}^{(\theta),\rm dist}={}&
-\frac12\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho\,
T_{\bar m\bar m,m}^{H,\theta}\right)
\right],
\\[1ex]
\widehat T_{m,4}^{(\theta),\rm dist}={}&
-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\bar\rho^2\Delta^{-1}
\widehat{\mathcal L}_{-1}^{[m,\theta]}
\left(\rho^{-2}\bar\rho^{-2}
T_{\bar mn,m}^{H,\theta}\right)
\right].
\end{aligned}
$$

定义

$$
S_{m,A}^{(\theta),\rm dist}
=-\frac{16\pi\Sigma}{\Delta^2R}
\widehat T_{m,A}^{(\theta),\rm dist}.
$$

Fourier projection 后、angular-coordinate change 前的完整中间方程为

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,\theta)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{(\theta),\rm dist}
}.
$$

## D.3 第二步：$y=-\cos\theta$

### D.3.1 Coordinate、measure 与 delta rules

定义

$$
\boxed{
y=-\cos\theta,
\qquad
D(y)=1-y^2,
\qquad
\sin\theta=\sqrt D,
\qquad
\partial_\theta=\sqrt D\,\partial_y
}.
$$

North-axis endpoint 是 $y=-1$，south-axis endpoint 是 $y=+1$。在 open chart $-1<y<1$，

$$
\boxed{
\sqrt{-g_{TRy\Phi}}
=\sqrt{-g_{TR\theta\Phi}}
\left|\frac{d\theta}{dy}\right|
=\frac{\Sigma r^2}{L^2}
},
$$

$$
\boxed{
\delta(\theta-\theta_p)
=\sin\theta_p\,\delta(y-y_p),
\qquad
y_p=-\cos\theta_p
}.
$$

因此 $\mathcal W_m^{(\theta)}$ 中的 $1/\sin\theta_p$ 与 angular delta Jacobian 精确相消，不能在后续再次加入。

### D.3.2 Vector、covector 与 Newman--Penrose rules

Angular components 按

$$
\boxed{
u^y=\sin\theta\,u^\theta,
\qquad
v_y=\frac{v_\theta}{\sin\theta}
}
$$

变换，所以 $N_p=n_Au_p^A$ 与 $\bar M_p=\bar m_Au_p^A$ 作为 scalars 不变。背景量变为

$$
\boxed{
\Sigma=r^2+a^2y^2,
\qquad
\rho=(r+iay)^{-1},
\qquad
\bar\rho=(r-iay)^{-1}
}.
$$

Kinnersley covectors 在 $(T,R,y,\Phi)$ 中为

$$
\boxed{
n_A^{\rm K}
=\left(
-\frac{\Delta}{2\Sigma},
\frac{\Delta r^2}{L^2\Sigma}
\left(1+\frac{2M}{r}\right),
0,
\frac{a\Delta D}{2\Sigma}
\right)
},
$$

$$
\boxed{
\bar m_A^{\rm K}
=\frac1{\sqrt2(r+iay)}
\left(
ia\sqrt D,
-\frac{iar^2\sqrt D}{L^2}
\left(2+\frac{4M}{r}\right),
\frac{\Sigma}{\sqrt D},
-i(r^2+a^2)\sqrt D
\right)
}.
$$

第二式是 open-chart component formula；axis behavior 由 spin-weighted regularity 完成。

### D.3.3 $y$ operators

Field angular operator 变为

$$
\boxed{
\mathcal A_{s,m}[f]
=D f_{,yy}-2y f_{,y}
-\frac{(m-sy)^2}{D}f+sf
}.
$$

Source operators 是

$$
\boxed{
\widehat{\mathcal L}^{[m]}_s
=\sqrt D\,\partial_y
-ia\sqrt D\,\partial_T
+\frac{m-sy}{\sqrt D}
},
$$

$$
\boxed{
\widehat{\mathcal J}_+
=-\left(2+\frac{4M}{r}\right)\partial_T
-\frac{L^2}{r^2}\partial_R
}.
$$

<a id="a4-y-operator"></a>

## D.4 完整 $(T,R,y;m)$ field operator

对任意 $s$，

$$
\begin{aligned}
\mathcal O_{s,m}^{(T,R,y)}[\psi_m]={}&
\left[
8M\left(2M-\frac{a^2R}{L^2}\right)
\left(1+\frac{2MR}{L^2}\right)-a^2D
\right]\partial_T^2\psi_m
\\
&-2\left[
L^2-(8M^2-a^2)\frac{R^2}{L^2}
+\frac{4a^2MR^3}{L^4}
\right]\partial_T\partial_R\psi_m
\\
&-\left(L^2-2MR+\frac{a^2R^2}{L^2}\right)
\frac{R^2}{L^2}\partial_R^2\psi_m
-\mathcal A_{s,m}[\psi_m]
\\
&+2iam\left(1+\frac{4MR}{L^2}\right)\partial_T\psi_m
+\frac{2iamR^2}{L^2}\partial_R\psi_m
\\
&+2\left[
2M\left(-s+2(s+2)\frac{MR}{L^2}
-3\frac{a^2R^2}{L^4}\right)
-\frac{a^2R}{L^2}-isay
\right]\partial_T\psi_m
\\
&+2R\left[-(1+s)+(s+3)\frac{MR}{L^2}
-2\frac{a^2R^2}{L^4}\right]\partial_R\psi_m
\\
&+\frac{2iamR}{L^2}\psi_m
+2\left[(1+s)\frac{MR}{L^2}
-\frac{a^2R^2}{L^4}\right]\psi_m.
\end{aligned}
$$

W02 的完整左端是

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,y)}
\equiv\left.\mathcal O_{s,m}^{(T,R,y)}\right|_{s=-2}
}.
$$

没有对 $\psi_{4,m}$ 作额外 angular rescaling。

<a id="a4-y-source"></a>

## D.5 三个 $y$-mode projections 与四块 exact source

### D.5.1 Mode projections

定义

$$
\boxed{
\mathcal W_m^{(y)}
=\frac{\mu L^2}
{2\pi\Sigma_pr_p^2u_p^T}
e^{-im\Phi_p(T)}
\delta(R-R_p(T))\delta(y-y_p(T))
}.
$$

于是

$$
\boxed{
T_{nn,m}^{H,y}=\mathcal W_m^{(y)}N_p^2,
\qquad
T_{\bar mn,m}^{H,y}=\mathcal W_m^{(y)}\bar M_pN_p,
\qquad
T_{\bar m\bar m,m}^{H,y}=\mathcal W_m^{(y)}\bar M_p^2
}.
$$

对 W02 赤道 orbit，$y_p=0$、$u_p^y=0$、$\Sigma_p=r_p^2$；$N_p,\bar M_p$ 与 Section V/附录 C 的 scalar contractions 相同。

### D.5.2 四个 ordered blocks

$$
\begin{aligned}
\widehat T_{m,1}^{\rm dist}={}&
-\rho^8\bar\rho\,
\widehat{\mathcal L}^{[m]}_{-1}
\left[
\rho^{-4}\widehat{\mathcal L}^{[m]}_0
\left(\rho^{-2}\bar\rho^{-1}T_{nn,m}^{H,y}\right)
\right],
\\[1ex]
\widehat T_{m,2}^{\rm dist}={}&
-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal L}^{[m]}_{-1}
\left[
\rho^{-4}\bar\rho^2\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho^{-2}\Delta^{-1}
T_{\bar mn,m}^{H,y}\right)
\right],
\\[1ex]
\widehat T_{m,3}^{\rm dist}={}&
-\frac12\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho\,
T_{\bar m\bar m,m}^{H,y}\right)
\right],
\\[1ex]
\widehat T_{m,4}^{\rm dist}={}&
-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\bar\rho^2\Delta^{-1}
\widehat{\mathcal L}^{[m]}_{-1}
\left(\rho^{-2}\bar\rho^{-2}
T_{\bar mn,m}^{H,y}\right)
\right].
\end{aligned}
$$

Block 2 是 outer $\widehat{\mathcal L}$ / inner $\widehat{\mathcal J}$；block 4 是 outer $\widehat{\mathcal J}$ / inner $\widehat{\mathcal L}$。不能交换。

定义

$$
\boxed{
S_{m,A}^{\rm dist}
=-\frac{16\pi\Sigma}{\Delta^2R}
\widehat T_{m,A}^{\rm dist}
}.
$$

最终 exact pre-Gaussian equation 是

$$
\boxed{
\mathcal O_{-2,m}^{(T,R,y)}[\psi_{4,m}]
=\sum_{A=1}^4S_{m,A}^{\rm dist}
}.
$$

外部 $\Sigma/(\Delta^2R)$ 是 field-point multiplier；$N_p,\bar M_p,r_p,\Phi_p$ 是 worldline functions；每个 nested operator 作用于其右侧括号中的全部 background factors、worldline amplitudes 和 exact distributions。

<a id="a4-weak-poles"></a>

## D.6 Dirac source 的分布积分意义

对 compactly supported smooth test function $\chi(T,R,y)$ 和

$$
F=C(T)\delta(R-R_p(T))\delta(y-y_p(T)),
$$

定义

$$
\langle F,\chi\rangle
=\int C(T)\chi(T,R_p(T),y_p(T))\,dT.
$$

这里的 derivative 是 distributional derivative。Coordinate derivatives 满足

$$
\boxed{
\langle\partial_TF,\chi\rangle
=-\langle F,\partial_T\chi\rangle,
\quad
\langle\partial_RF,\chi\rangle
=-\langle F,\partial_R\chi\rangle,
\quad
\langle\partial_yF,\chi\rangle
=-\langle F,\partial_y\chi\rangle
}.
$$

一般一阶 operator

$$
\mathcal D=a^T\partial_T+a^R\partial_R+a^y\partial_y+b
$$

的 bilinear transpose 是

$$
\mathcal D^{\rm tr}\chi
=-\partial_T(a^T\chi)-\partial_R(a^R\chi)
-\partial_y(a^y\chi)+b\chi.
$$

这是 integration by parts 的定义：每把一个 derivative 从 distribution 移到 test function 都产生一个负号。Nested operators 的 transposes 必须以反序作用，所以该定义不仅固定 delta-derivative signs，也固定两个 mixed source blocks 的 ordering。它自动包含 moving support 与 $e^{-im\Phi_p(T)}$ 的时间依赖；本阶段不把它们展开成 explicit delta jets。

## D.7 旋转轴 $y=\pm1$ 的 spin-weighted regularity

这里的“轴端”指几何旋转轴

$$
y=-1\Longleftrightarrow\theta=0,
\qquad
y=+1\Longleftrightarrow\theta=\pi,
$$

不是复平面中的极点。$\mathcal A_{s,m}$ 含 $(1-y^2)^{-1}$，但这只是使用 azimuthal spin frame 后的表面 coordinate singularity；regular spin-weighted mode 必须以相应阶数在轴上消失。由 singular part 的 indicial equation，north 与 south 的 regular powers 是

$$
\boxed{
p_- =\frac{|m+s|}{2},
\qquad
p_+ =\frac{|m-s|}{2}
}.
$$

因此

$$
\psi_{4,m}
\sim(1+y)^{p_-}(1-y)^{p_+}
\times(\text{axis-smooth function}),
$$

且对 $s=-2$：

| $m$ | $p_-$ at $y=-1$ | $p_+$ at $y=+1$ | $\theta$ powers |
|---:|---:|---:|---|
| 2 | 0 | 2 | $\theta^0$ north，$(\pi-\theta)^4$ south |
| 4 | 1 | 3 | $\theta^2$ north，$(\pi-\theta)^6$ south |

这些是 continuum field regularity conditions，不是对 finite-width angular extension 的定义，也尚未给出数值 stencil/ghost/parity closure。显式 Wigner-$d$ functions 与这些 powers 一致，但该检查只作为 axis regularity 的辅助证据；A4 source transformation 的主结论是 Fourier/delta/Jacobian/operators/four-block endpoint 的一致性。

## D.8 Evidence、等价性与边界

既有 [A4 result](../../../../data/kerr_point_particle_verification/results/A4_fixed_m_y_distribution.json) 已给出：

- $\partial_\Phi\to im$、field angular operator、$\widehat{\mathcal L}_s^{[m]}$、$\widehat{\mathcal J}_+$ residuals 为 $0$；
- 完整 $\theta\to y$ field-operator comparison residual 为 $0$；
- 把三个 Kinnersley projections 视为 arbitrary functions 时，四个 ordered source blocks 和 total residual 为 $0$；
- Fourier/delta/weight/Jacobian、vector/covector contraction 与分布求导符号分别检查通过；
- Wigner-$d$ axis powers 与上述 $(0,2)$、$(1,3)$ 一致；这是辅助 regularity check，不替代 numerical angular boundary closure。

升级后的 staged A4 从 Section V 的完整未投影 $E_0$ 出发。生成路线先让未投影 operators 作用于 $e^{im(\Phi-\Phi_p)}/(2\pi)$ 和完整点粒子 operands，再抽取 $e^{im\Phi}$；理论 $\theta$/$y$ endpoints 只出现在比较侧。结果为

$$
\operatorname{Res}E_m^{(\theta)}=0,
\qquad
\operatorname{Res}S_{m,A}^{(\theta),\rm dist}=0,
\quad A=1,\ldots,4,
$$

$$
\operatorname{Res}E_m^{(y)}=0,
\qquad
\operatorname{Res}S_{m,A}^{\rm dist}=0,
\quad A=1,\ldots,4.
$$

$y$ 路线还分别比较了 angular-delta Jacobian、$\Sigma(R,\theta)\to\Sigma(R,y)$ 外因子和四个完整 source blocks；生成结果不含残余 $\theta$ 或 $\Phi$。较大的源差式在独立 worldline/distributional-jet 基上展开，最大约 $8.7\times10^5$ leaves，全部严格为 $0$。可重放代码和结果见 [A4 checks](../../../../code/kerr_point_particle_verification/checks/a4/)与 [A4 results](../../../../data/kerr_point_particle_verification/results/a4/)；汇总给出 `A4FullEndpointVerified=True`。旧 arbitrary-projection checks 与轴正则性检查仅作为辅助证据，不用于关闭完整端点。

开放边界：

- A1 signature/Weyl-scalar、sourced-unit normalization 与 independent Ripley-PDF transcription 仍开放；
- A2/A3 endpoints 已验证，A4 继承其 assumptions；
- A4 完整 $\theta$ stage 与 full point-particle assembly 已验证；
- interior identities 使用 $R>0$、$-1<y<1$、$\Delta\ne0$，particle 不位于 angular coordinate axis endpoint；
- $y=\pm1$ 的 continuum spin-weighted powers 已识别，但 numerical angular boundary condition 尚未闭合；
- 不包含 Gaussian、off-worldline extension、time jets、boundary limits 或 discretization。
