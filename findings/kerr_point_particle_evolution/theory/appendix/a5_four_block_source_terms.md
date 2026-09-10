# A5 四块 Gaussian 点粒子数值源

- 主文：[DERIVATION.md Section VII](../DERIVATION.md#vii-gaussianradial-extension-与-time-jets)
- 完整 jets：[附录 E](W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md)
- 上游 exact source：[附录 D.5](W02_1_4_FIXED_M_Y_DERIVATION.md#a4-y-source)
- 验证证据：[A5 endpoint summary](../../../../data/kerr_point_particle_verification/results/a5/a5_endpoint_summary.json)
- 状态：`A5FullEndpointVerifiedWithinInteriorScope=True`
- 范围：fixed widths、$0<R<R_H$、$|y|<1$ 内的 time algebra；全部 coordinate $\partial_T$ 已消去，$\partial_R,\partial_y$ 保留给空间离散
- 开放边界：finite-domain Gaussian normalization、$y=\pm1$ axis completion、$R=0,R_H$ endpoint rows、turning branch 与 finite-width convergence

## 共同定义

### Background 与空间算符系数

背景参数为 Kerr $(M,a)$、compactification scale $L$、particle mass $\mu$ 和 fixed azimuthal mode $m$。令

$$
r=\frac{L^2}{R},
\qquad
\widehat\Delta(R)=R^2\Delta
=L^4-2ML^2R+a^2R^2,
$$

$$
\widehat\Sigma(R,y)=R^2\Sigma
=L^4+a^2R^2y^2,
$$

$$
\Delta=\frac{\widehat\Delta}{R^2},
\qquad
\Sigma=\frac{\widehat\Sigma}{R^2},
\qquad
r_+=M+\sqrt{M^2-a^2},
\qquad
R_H=\frac{L^2}{r_+},
$$

$$
\rho(R,y)=\frac{R}{L^2+iaRy},
\qquad
\bar\rho(R,y)=\frac{R}{L^2-iaRy},
$$

$$
c(y)=\sqrt{1-y^2},
\qquad
\ell_s(y)=\frac{m-sy}{c(y)},
$$

其中 $s=0$ 只用于 block 1 的 inner $L_0$；blocks 1、2 的 outer $L_{-1}$ 和 block 4 的 inner $L_{-1}$ 都使用 $s=-1$。

$$
j_T(R)=-\left(2+\frac{4MR}{L^2}\right),
\qquad
j_R(R)=-\frac{R^2}{L^2}.
$$

本文所有 products 都是 field-point multiplication；$\partial_R$、$\partial_y$ 作用于其右侧方括号内的完整 operand。

### Fixed-width Gaussian

取固定 $\sigma_R,\sigma_y>0$，并定义

$$
G_R(T,R)=\frac1{\sqrt{2\pi}\sigma_R}
\exp\left[-\frac{(R-R_p(T))^2}{2\sigma_R^2}\right],
$$

$$
G_y(T,y)=\frac1{\sqrt{2\pi}\sigma_y}
\exp\left[-\frac{(y-y_p(T))^2}{2\sigma_y^2}\right],
\qquad
G=G_RG_y.
$$

它们相对于 full-line coordinate measures $dR,dy$ 单位归一。W02 的赤道 orbit 满足

$$
y_p=0,
\qquad
\dot y_p=\ddot y_p=0.
$$

记

$$
x=R-R_p(T),
$$

$$
\Gamma_1=\frac{x\dot R_p}{\sigma_R^2},
$$

$$
\Gamma_2
=\frac{x^2\dot R_p^2}{\sigma_R^4}
+\frac{x\ddot R_p-\dot R_p^2}{\sigma_R^2}.
$$

### Amplitude、phase 与 Gaussian-dressed jets

使用 horizon-regular amplitudes

$$
\mathcal A_{nn}=\widetilde A_{nn}=A_0\widehat N_p^2,
$$

$$
\mathcal A_{\bar mn}=\widetilde A_{\bar mn}
=A_0\bar M_p\widehat N_p,
$$

$$
\mathcal A_{\bar m\bar m}=A_{\bar m\bar m}
=A_0\bar M_p^2,
$$

其中

$$
A_0=\frac{\mu L^2}{2\pi r_p^4u_p^T}.
$$

$R_p,\dot R_p,\ddot R_p,\Phi_p,\dot\Phi_p,\ddot\Phi_p$，以及 $A_0,\widehat N_p,\bar M_p$ 和上述三个 amplitudes 的完整 $0,1,2$ 阶 jets，均由 [附录 E.4--E.5](W02_1_5_GAUSSIAN_EXTENSION_TIME_JETS.md#a5-worldline-jets) 定义并已在 A5 验证；它们不是本页的额外输入。

对 $ab\in\{nn,\bar mn,\bar m\bar m\}$，定义

$$
Q_{ab}^{(0)}=\mathcal A_{ab}e^{-im\Phi_p},
$$

$$
Q_{ab}^{(1)}
=e^{-im\Phi_p}
\left(\dot{\mathcal A}_{ab}
-im\dot\Phi_p\mathcal A_{ab}\right),
$$

$$
Q_{ab}^{(2)}
=e^{-im\Phi_p}
\left[
\ddot{\mathcal A}_{ab}
-2im\dot\Phi_p\dot{\mathcal A}_{ab}
-im\ddot\Phi_p\mathcal A_{ab}
-m^2\dot\Phi_p^2\mathcal A_{ab}
\right].
$$

三个 Gaussian-dressed jet sequences 统一定义为

$$
f_{ab}^{(0)}=Q_{ab}^{(0)}G,
$$

$$
f_{ab}^{(1)}
=\left(Q_{ab}^{(1)}+Q_{ab}^{(0)}\Gamma_1\right)G,
$$

$$
f_{ab}^{(2)}
=\left(Q_{ab}^{(2)}
+2Q_{ab}^{(1)}\Gamma_1
+Q_{ab}^{(0)}\Gamma_2\right)G.
$$

与附录 E 的 notation 对应关系是：本页 $f_{nn}^{(k)}$ 即 $\widetilde f_{nn}^{(k)}$，$f_{\bar mn}^{(k)}$ 即 $\widetilde f_{\bar mn}^{(k)}$；$\bar m\bar m$ sequence 不带 tilde，因为该 projection 不含 compactified-$\Delta$ rescaling。下文统一使用本页的 $f_{ab}^{(k)}$。

Off-worldline projections 是

$$
F_{nn}^{\rm ext}=\widehat\Delta^2 f_{nn}^{(0)},
\qquad
F_{\bar mn}^{\rm ext}=\widehat\Delta f_{\bar mn}^{(0)},
\qquad
F_{\bar m\bar m}^{\rm ext}=f_{\bar m\bar m}^{(0)}.
$$

下面的 formulas 已把这些 $\widehat\Delta$ factors 与 A4 blocks 中的 $\Delta^{-1},\Delta^2$ 及外部 $-16\pi\Sigma/(\Delta^2R)$ 逐块解析组合：block 1 的 $\widehat\Delta^2$ 穿过两个 $L$ factors 后与外部 $\widehat\Delta^{-2}$ 相消；block 2 的 inner $\Delta^{-1}\widehat\Delta$ 先成为 $R^2$，剩余 outer $\Delta^2$ 再与外部因子相消；block 3 的 outer $\Delta^2$ 直接与外部因子相消；block 4 的 $\widehat\Delta$ 穿过 inner $L_{-1}$ 后与 $\Delta^{-1}$ 形成 $R^2$，outer $\Delta^2$ 再相消。不得在实现中重新乘入这些 factors。

## Block 1：$nn$，LL ordering

对 $k=0,1,2$，

$$
U_1^{(k)}
=\rho^{-2}\bar\rho^{-1}f_{nn}^{(k)}.
$$

Inner $L_0$ 给出

$$
V_1^{(0)}
=c\,\partial_yU_1^{(0)}
-iac\,U_1^{(1)}+\ell_0U_1^{(0)},
$$

$$
V_1^{(1)}
=c\,\partial_yU_1^{(1)}
-iac\,U_1^{(2)}+\ell_0U_1^{(1)}.
$$

定义

$$
W_1^{(0)}=\rho^{-4}V_1^{(0)},
\qquad
W_1^{(1)}=\rho^{-4}V_1^{(1)}.
$$

Outer $L_{-1}$ 后

$$
\boxed{
\mathfrak H_1
=-\rho^8\bar\rho
\left[
c\,\partial_yW_1^{(0)}
-iac\,W_1^{(1)}+\ell_{-1}W_1^{(0)}
\right]
},
$$

$$
\boxed{
S_{m,1}^{\rm num}
=-16\pi R\widehat\Sigma\,\mathfrak H_1
}.
$$

## Block 2：$\bar mn$，LJ ordering

Inner $J$ 前的 $\Delta^{-1}\widehat\Delta$ 已解析相消为 $R^2$。对 $k=0,1,2$，

$$
U_2^{(k)}
=\rho^{-2}\bar\rho^{-2}R^2f_{\bar mn}^{(k)}.
$$

Inner $J$ 给出

$$
V_2^{(0)}
=j_TU_2^{(1)}+j_R\partial_RU_2^{(0)},
$$

$$
V_2^{(1)}
=j_TU_2^{(2)}+j_R\partial_RU_2^{(1)}.
$$

定义

$$
W_2^{(0)}=\rho^{-4}\bar\rho^2V_2^{(0)},
\qquad
W_2^{(1)}=\rho^{-4}\bar\rho^2V_2^{(1)}.
$$

Outer $L_{-1}$ 后

$$
\boxed{
\mathfrak H_2
=-\frac{\rho^8\bar\rho}{\sqrt2R^4}
\left[
c\,\partial_yW_2^{(0)}
-iac\,W_2^{(1)}+\ell_{-1}W_2^{(0)}
\right]
},
$$

$$
\boxed{
S_{m,2}^{\rm num}
=-16\pi R\widehat\Sigma\,\mathfrak H_2
}.
$$

## Block 3：$\bar m\bar m$，JJ ordering

对 $k=0,1,2$，

$$
U_3^{(k)}
=\rho^{-2}\bar\rho f_{\bar m\bar m}^{(k)}.
$$

Inner $J$ 给出

$$
V_3^{(0)}
=j_TU_3^{(1)}+j_R\partial_RU_3^{(0)},
$$

$$
V_3^{(1)}
=j_TU_3^{(2)}+j_R\partial_RU_3^{(1)}.
$$

定义

$$
W_3^{(0)}=\rho^{-4}V_3^{(0)},
\qquad
W_3^{(1)}=\rho^{-4}V_3^{(1)}.
$$

Outer $J$ 后

$$
\boxed{
\mathfrak H_3
=-\frac{\rho^8\bar\rho}{2R^4}
\left[
j_TW_3^{(1)}+j_R\partial_RW_3^{(0)}
\right]
},
$$

$$
\boxed{
S_{m,3}^{\rm num}
=-16\pi R\widehat\Sigma\,\mathfrak H_3
}.
$$

## Block 4：$\bar mn$，JL ordering

对 $k=0,1,2$，

$$
U_4^{(k)}
=\rho^{-2}\bar\rho^{-2}f_{\bar mn}^{(k)}.
$$

Inner $L_{-1}$ 给出

$$
V_4^{(0)}
=c\,\partial_yU_4^{(0)}
-iac\,U_4^{(1)}+\ell_{-1}U_4^{(0)},
$$

$$
V_4^{(1)}
=c\,\partial_yU_4^{(1)}
-iac\,U_4^{(2)}+\ell_{-1}U_4^{(1)}.
$$

Inner $L$ 让 field-point-only $\widehat\Delta$ 穿过，随后与 $\Delta^{-1}$ 相消为 $R^2$。定义

$$
W_4^{(0)}=\rho^{-4}\bar\rho^2R^2V_4^{(0)},
\qquad
W_4^{(1)}=\rho^{-4}\bar\rho^2R^2V_4^{(1)}.
$$

Outer $J$ 后

$$
\boxed{
\mathfrak H_4
=-\frac{\rho^8\bar\rho}{\sqrt2R^4}
\left[
j_TW_4^{(1)}+j_R\partial_RW_4^{(0)}
\right]
},
$$

$$
\boxed{
S_{m,4}^{\rm num}
=-16\pi R\widehat\Sigma\,\mathfrak H_4
}.
$$

## 总源与数组映射

$$
\boxed{
S_m^{\rm num}
=S_{m,1}^{\rm num}+S_{m,2}^{\rm num}
+S_{m,3}^{\rm num}+S_{m,4}^{\rm num}
}.
$$

对于后续 Python finite-difference implementation，可把每个 $U_A^{(k)},V_A^{(k)},W_A^{(k)},\mathfrak H_A,S_{m,A}^{\rm num}$ 表示为同一 $(n_R,n_y)$ complex array。Field-point multiplication 是逐元素乘法；$D_R$ 只替代显式 $\partial_R$，$D_y$ 只替代显式 $\partial_y$；每块必须由 inner 到 outer 计算。上述 formulas 是 interior expressions，不能在 $R=0,R_H$ 或 $y=\pm1$ rows 直接依赖 $R^{-4}$、$c^{-1}$ 的浮点抵消；这些 rows 需要单独的数值边界 closure。
