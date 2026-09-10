# 附录 A：$(T,R,\theta,\Phi)$ 中的 Ripley 目标算符

- Parent: [DERIVATION.md](../DERIVATION.md)
- Scope: W02.1.1
- Status: conditional 12-jet algebra passed；independent PDF transcription remains open
- Primary source: Ripley et al., Eqs. (21)--(22)

## Evolved field

Ripley defines the regular-tetrad peeling field by

$$
\Psi_4^{\rm R}=R\psi_4.
$$

The equation below governs $\psi_4(T,R,\theta,\Phi)$. No azimuthal or $y=-\cos\theta$ projection has been made.

## Spin-weighted angular operator

Define

$$
{}_s\!\Delta_{S^2}f
=\frac1{\sin\theta}\partial_\theta
\left(\sin\theta\,\partial_\theta f\right)
+\left[
s-\frac{(-i\partial_\Phi+s\cos\theta)^2}{\sin^2\theta}
\right]f.
$$

## Ripley Eq. (22)

For spin weight $s$,

$$
\begin{aligned}
\mathcal O_{\rm R}[\psi_4]={}&
\left[
8M\left(2M-\frac{a^2R}{L^2}\right)
\left(1+\frac{2MR}{L^2}\right)
-a^2\sin^2\theta
\right]\partial_T^2\psi_4
\\
&-2\left[
L^2-(8M^2-a^2)\frac{R^2}{L^2}
+\frac{4a^2MR^3}{L^4}
\right]\partial_T\partial_R\psi_4
\\
&-\left(
L^2-2MR+\frac{a^2R^2}{L^2}
\right)\frac{R^2}{L^2}\partial_R^2\psi_4
-{}_s\!\Delta_{S^2}\psi_4
\\
&+2a\left(1+\frac{4MR}{L^2}\right)
\partial_T\partial_\Phi\psi_4
+\frac{2aR^2}{L^2}\partial_R\partial_\Phi\psi_4
\\
&+2\left[
2M\left(
-s+2(s+2)\frac{MR}{L^2}
-3\frac{a^2R^2}{L^4}
\right)
-\frac{a^2R}{L^2}
+isa\cos\theta
\right]\partial_T\psi_4
\\
&+2R\left[
-(1+s)+(s+3)\frac{MR}{L^2}
-2\frac{a^2R^2}{L^4}
\right]\partial_R\psi_4
\\
&+\frac{2aR}{L^2}\partial_\Phi\psi_4
+2\left[
(1+s)\frac{MR}{L^2}
-\frac{a^2R^2}{L^4}
\right]\psi_4.
\end{aligned}
$$

The vacuum target equation is

$$
\boxed{\mathcal O_{\rm R}[\psi_4]=0}.
$$

W02 uses $s=-2$. The compatible source is constructed in W02.1.2; this appendix does not yet set $\partial_\Phi\to im$.

## Scope and checks

This transcription was cross-read against the target repository's current coefficient implementation, but code agreement is not an independent proof. W02.2 must:

1. re-enter Eq. (22) from the PDF;
2. compare every zero-, first- and second-derivative coefficient;
3. verify the angular-operator sign independently;
4. only then derive the fixed-$m$, $y$ form in W02.1.4/W02.2.

The conditional 12-jet comparison has passed; see [A1/A2 code](../../../../code/kerr_point_particle_verification/checks/A1_A2_operator.wl) and [result](../../../../data/kerr_point_particle_verification/results/A1_A2_operator.json). Item 1 above and the A1 convention bridge remain open, so this is not yet an approved operator.
