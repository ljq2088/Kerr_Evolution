# 附录 B：Sasaki--Tagoshi source 算符与 hyperboloidal 拉回

- 主文：[DERIVATION.md](../DERIVATION.md)
- 范围：W02.1.2
- 状态：four-block pullback algebra 条件通过；A1 convention/unit bridge 仍开放
- 主要来源：Sasaki--Tagoshi Eqs. (6)--(12)

## 在推导中的作用

[DERIVATION.md](../DERIVATION.md) 给出从 Sasaki 方程到 Ripley-compatible source 的主线，本附录保留 W02.1.2 的完整公式。当前内容的逻辑顺序是：从原始 BL 有源方程及 $B$ 对象的定义出发，展开四个有序 source blocks，再进行坐标拉回并写出完整 hyperboloidal source，最后给出它与 Ripley 算符之间的整体关系。

## 原始 BL 方程

在 Sasaki--Tagoshi convention 下，

$$
\mathcal O_{\rm ST}[\phi_{\rm ST}]
=4\pi\Sigma\widehat T_{\rm ST},
\qquad
\phi_{\rm ST}=\rho^{-4}\Psi_4^{\rm K},
$$

$$
\widehat T_{\rm ST}=2(B'_2+B_2^{*\prime}).
$$

$T_{\mu\nu}$ 是 matter stress tensor。三个量

$$
T_{nn},\qquad T_{\bar mn},\qquad T_{\bar m\bar m}
$$

是它在 Kinnersley tetrad 上的 projections。$B'_2$ 和 $B_2^{*\prime}$ 是这些 projections 的微分组合，不是额外的 matter variables。

## 四个有序 source blocks

$$
\begin{aligned}
B'_2={}&
-\frac12\rho^8\bar\rho\,
\mathcal L_{-1}
\left[
\rho^{-4}\mathcal L_0
\left(\rho^{-2}\bar\rho^{-1}T_{nn}\right)
\right]\\
&-\frac{1}{2\sqrt2}\rho^8\bar\rho\,\Delta^2
\mathcal L_{-1}
\left[
\rho^{-4}\bar\rho^2\mathcal J_+
\left(
\rho^{-2}\bar\rho^{-2}\Delta^{-1}T_{\bar mn}
\right)
\right],
\end{aligned}
$$

$$
\begin{aligned}
B_2^{*\prime}={}&
-\frac14\rho^8\bar\rho\,\Delta^2
\mathcal J_+
\left[
\rho^{-4}\mathcal J_+
\left(\rho^{-2}\bar\rho T_{\bar m\bar m}\right)
\right]\\
&-\frac{1}{2\sqrt2}\rho^8\bar\rho\,\Delta^2
\mathcal J_+
\left[
\rho^{-4}\bar\rho^2\Delta^{-1}\mathcal L_{-1}
\left(
\rho^{-2}\bar\rho^{-2}T_{\bar mn}
\right)
\right].
\end{aligned}
$$

其中算符为

$$
\mathcal L_s
=\partial_\theta-i\csc\theta\,\partial_\phi
-ia\sin\theta\,\partial_t+s\cot\theta,
$$

$$
\mathcal J_+
=\partial_r
-\frac{r^2+a^2}{\Delta}\partial_t
-\frac a\Delta\partial_\phi.
$$

嵌套运算必须由内向外读取。特别地，两个混合的 $T_{\bar mn}$ blocks 具有不同的 operator ordering。$B_2^{*\prime}$ 是 Sasaki 对第二组项的命名，不能用对整个 $B'_2$ 作机械复共轭来替代。

依赖关系为

```text
T_{mu nu}
  -> {T_nn, T_mbar n, T_mbar mbar}
  -> four differential blocks
  -> {B'_2, B_2^{*prime}}
  -> hat T_ST
  -> right-hand side 4 pi Sigma hat T_ST
```

对于点粒子，delta distributions 通过三个 projections 进入；随后微分 blocks 产生 delta derivatives。

## 坐标拉回

采用

$$
T=t+h(r),\qquad R=\frac{L^2}{r},\qquad
\Phi=\phi+q(r),
$$

$$
h'=\frac{r^2+a^2}{\Delta}-2-\frac{4M}{r},
\qquad q'=\frac a\Delta.
$$

derivative map 为

$$
\partial_t=\partial_T,
\qquad
\partial_\phi=\partial_\Phi,
$$

$$
\partial_r\big|_{t,\phi}
=h'\partial_T-\frac{L^2}{r^2}\partial_R
+\frac a\Delta\partial_\Phi.
$$

因此

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

在进行任何 fixed-$m$ projection 之前，$\widehat{\mathcal J}_+$ 中的 $\partial_\Phi$ 项已经解析相消。

## 拉回后的 source

把外部因子 2 乘入四个 blocks 后，

$$
\begin{aligned}
\widehat T_H={}&
-\rho^8\bar\rho\,
\widehat{\mathcal L}_{-1}
\left[
\rho^{-4}\widehat{\mathcal L}_0
\left(\rho^{-2}\bar\rho^{-1}T_{nn}^H\right)
\right]\\
&-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal L}_{-1}
\left[
\rho^{-4}\bar\rho^2\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho^{-2}\Delta^{-1}T_{\bar mn}^H\right)
\right]\\
&-\frac12\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\widehat{\mathcal J}_+
\left(\rho^{-2}\bar\rho T_{\bar m\bar m}^H\right)
\right]\\
&-\frac1{\sqrt2}\rho^8\bar\rho\,\Delta^2
\widehat{\mathcal J}_+
\left[
\rho^{-4}\bar\rho^2\Delta^{-1}
\widehat{\mathcal L}_{-1}
\left(\rho^{-2}\bar\rho^{-2}T_{\bar mn}^H\right)
\right].
\end{aligned}
$$

所有光滑因子均视为 $(T,R,\theta,\Phi)$ 的函数，并取 $r=L^2/R$。微分算符同时作用于光滑因子和 distributions。

## 与 Ripley 算符的候选关系

目标仓库给出的关系为

$$
\mathcal O_{\rm ST}
\left[\frac{\Delta^2R}{4}\psi_4\right]
=-\frac{\Delta^2R}{4}\mathcal O_{\rm R}[\psi_4].
$$

在该 convention bridge 下，

$$
-\frac{\Delta^2R}{4}\mathcal O_{\rm R}[\psi_4]
=4\pi\Sigma\widehat T_H,
$$

因此

$$
\boxed{
\mathcal O_{\rm R}[\psi_4]
=-\frac{16\pi\Sigma}{\Delta^2R}\widehat T_H
}.
$$

状态：当前仓库中的 A2 operator/four-block residuals 均为零，见 [operator 结果](../../../../data/kerr_point_particle_verification/results/A1_A2_operator.json)和 [four-block 结果](../../../../data/kerr_point_particle_verification/results/A2_four_source_blocks.json)。signature/unit convention bridge 仍开放，因此该有源方程仍为条件成立。

## 范围边界

本附录尚未进行 fixed-$m$、$y=-\cos\theta$、Gaussian regularization 或 source time-derivative elimination；这些工作从 W02.1.4 开始。
