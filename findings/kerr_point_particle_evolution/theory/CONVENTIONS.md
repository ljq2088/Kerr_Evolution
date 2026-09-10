# W02.1 conventions

- Scope: W02.1.1--W02.1.3
- Status: A2/A3 algebra conditionally passed；A1 convention gate remains open
- Last updated: 2026-08-31

## Purpose

This file is the single notation ledger for the current theory candidate. It distinguishes the Sasaki--Tagoshi source convention from the Ripley evolution convention instead of silently identifying them.

## Evidence labels

- `primary`: transcribed from the cited paper.
- `derived`: obtained by an explicit analytic transformation recorded in `DERIVATION.md`.
- `legacy_verified`: the target repository records a Mathematica residual, but it has not been replayed in this repository.
- `awaiting_verification`: requires W02.2 before use as an approved equation.

## Background and units

Use $c=1$. For the sourced candidate, retain Sasaki's displayed $4\pi\Sigma\widehat T$ normalization, corresponding to its gravitational-unit convention. Ripley states $8\pi G=1$, but the equation borrowed from Ripley is vacuum; the sourced unit bridge is therefore not inferred from Ripley's vacuum equation and remains a W02.2 normalization check.

The derivation retains

$$
M,\qquad a,\qquad L,\qquad \mu,
$$

with $|a|<M$ in the nonextremal domain considered here. Define

$$
\Sigma=r^2+a^2\cos^2\theta,
\qquad
\Delta=r^2-2Mr+a^2=(r-r_+)(r-r_-),
$$

$$
r_\pm=M\pm\sqrt{M^2-a^2}.
$$

$\rho=(r-ia\cos\theta)^{-1}$ is the Newman--Penrose quantity, not a matter density.

## Two paper conventions

### Sasaki--Tagoshi source convention

Sasaki--Tagoshi Eqs. (3)--(12) use the Kerr line element with signature $(-,+,+,+)$ and define

$$
\Psi_4^{\rm K}
=-C_{\alpha\beta\gamma\delta}
n_{\rm K}^\alpha\bar m_{\rm K}^\beta
n_{\rm K}^\gamma\bar m_{\rm K}^\delta,
$$

$$
\phi_{\rm ST}=\rho^{-4}\Psi_4^{\rm K}.
$$

The source is built from Kinnersley projections

$$
T_{nn},\qquad T_{\bar mn},\qquad T_{\bar m\bar m}.
$$

Status: `primary`.

### Ripley evolution convention

Ripley et al. explicitly state signature $(+,-,-,-)$ and use a boost/spin-rotated tetrad regular at the future horizon and SCRI+. Their peeling variable is

$$
\Psi_4^{\rm R}=R\psi_4.
$$

The target evolution operator is Ripley Eq. (22), transcribed in unprojected $(T,R,\theta,\Phi)$ form in [the operator appendix](appendix/RIPLEY_OPERATOR.md).

Status: `primary` transcription, awaiting independent W02.2 coefficient check.

### Signature bridge

The target repository uses

$$
\rho^{-4}\Psi_4^{\rm K}
=\frac{\Delta^2}{4}\Psi_4^{\rm R}
$$

after applying the Ripley tetrad boost/spin rotation.

In the current repository `A1MasterFieldResidual = 0` conditionally verifies the entered boost/spin algebra; see [A1 code](../../../code/kerr_point_particle_verification/checks/A1_A2_operator.wl) and [result](../../../data/kerr_point_particle_verification/results/A1_A2_operator.json). The same execution explicitly leaves the signature/Weyl-scalar and sourced-unit bridges undecided. A sign change there would propagate to the source normalization.

## Coordinates

Use BL coordinates $(t,r,\theta,\phi)$ and hyperboloidal coordinates $(T,R,\theta,\Phi)$:

$$
T=t+h(r),\qquad R=\frac{L^2}{r},\qquad \Phi=\phi+q(r),
$$

$$
h(r)=r_*-2r-4M\ln(r/M),
$$

$$
\frac{dr_*}{dr}=\frac{r^2+a^2}{\Delta},
\qquad
q'(r)=\frac a\Delta.
$$

The additive constants in $r_*$ and $q$ are fixed once per run and used consistently for field and worldline data. $R=0$ is SCRI+ and $R_H=L^2/r_+$ is the future horizon.

With the chosen asymptotic constant,

$$
T\big|_{\mathscr I^+}=u=t-r_*.
$$

Status: `primary` for the Ripley coordinate choice; local derivative maps are `derived`.

## Fourier and angular conventions

The future fixed-$m$ decomposition will use

$$
\psi_4(T,R,\theta,\Phi)
=\sum_m\psi_{4,m}(T,R,\theta)e^{im\Phi}.
$$

Hence a point source carries $e^{-im\Phi_p}$. W02.1.1--1.3 do not yet perform this projection.

The numerical angular coordinate will be $y=-\cos\theta$, but that transformation belongs to W02.1.4 and is outside the present execution scope.

## Source notation

$T^{AB}_H$ is the matter stress tensor in hyperboloidal coordinates. $T_{nn}^H,T_{\bar mn}^H,T_{\bar m\bar m}^H$ are its Kinnersley projections. $B'_2$ and $B_2^{*\prime}$ are differential combinations of those projections. Finally,

$$
\widehat T_H=2(B'_2+B_2^{*\prime})\big|_{\rm pulled\ back}
$$

is the scalar source entering the transformed Teukolsky equation. These objects must not be denoted by the same unqualified symbol $T$.

## Current scope exclusions

Not yet defined here:

- fixed-$m$ and $y$ projections;
- Gaussian regularization and off-worldline extension;
- analytic source time jets;
- first-order reduction and numerical boundaries;
- strain conversion.

Those are not silently assumed in W02.1.1--1.3.
