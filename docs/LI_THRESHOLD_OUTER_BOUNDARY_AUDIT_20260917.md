# Near-threshold outer-boundary audit: fixed fresh source

The currently used Coulomb/Whittaker boundary is stable under `Rout=32000→64000M`: the scalar22 field changes by less than `3e-9` at the five prescribed radii. A deliberately premature large-radius series boundary can nevertheless move the bound-field nodes strongly, and can separately corrupt the conversion from a computed coefficient to an infinity flux. This establishes a concrete numerical failure class. It does **not** establish the original authors' outer cutoff or prove that this caused their plotted discrepancy.

All source samples are the fresh `r_p=41.1,42.1`, scalar22, `Lmetric=6`, `ntheta=12`, `nr=8`, horizon-log-order32 results. The source support, source values, source normalization, cloud frequency, and other 17 scalar channels are fixed. The background remains the stationary synchronized `a=.8771530275949366`, not a newly calculated `a=.88` environment. No production equation or solver file was changed. All outer radii were prescribed before evaluating agreement with reference pixels.

## What the paper actually specifies

The local Li v2 `source/main.tex` lines 592–604 explicitly prescribe

\[
P(r)=\frac{e^{ik r_*}r^{i\mu^2/k}}{\sqrt{r^2+a^2}},\qquad
R_{\rm up}=P(r)\sum_{j=0}^{4}B_jr^{-j},\quad B_0=1.
\]

The expansion order is therefore **known: four**. The downloaded text does not provide a numerical outer endpoint. Line 729 lists cutoff and boundary choices among possible numerical differences. The bibliography links the coefficient notebook at `https://github.com/dongjun826/EMRI-in-scalar-clouds.git`; the separate reference audit fixes its commit and files. The public artifact has coefficient expressions, not a complete integration driver that reveals `Rout`.

The first control uses leading-asymptotic and sixth-order series in `exp(ikr) r^(-1+i beta/k)` at `1000,4000,8000M`. These are sensitivity controls, not Li's actual finite-order prescription. A second control implements the **exact published prefactor** and independently derives its four coefficients at `4000,8000,32000M`.

For the second calculation,

\[
\frac{P'}P=ik\frac{r^2+a^2}{\Delta}+\frac{i\mu^2}{kr}-\frac{r}{r^2+a^2}.
\]

Substituting `R=P S`, `S=sum B_j x^j`, `x=1/r`, into the polynomial radial ODE gives a triangular coefficient recurrence; the coefficient multiplying the new `B_n` is `-2 i k n` at power `x^(n+1)`. This was evaluated with 65-digit arithmetic, independently of `environment_radial.infinity_series`. The reference audit independently evaluated the original-author notebook B1–B4 expressions for both the propagating and bound cases: all eight coefficients agree with this recurrence, with maximum relative difference `1.4333e-45`, limited by the retained 45-digit strings. Evidence is `docs/field_alignment_20260917/li_author_boundary_coefficients_check.json`. This is evaluation of the author’s exact coefficient expressions through a restricted mathematical parser and80-digit mpmath, **not execution of the original Mathematica notebook**. It establishes that the present fourth-order prefactor/coefficients correspond to the published ones; it does not reveal the author’s numerical outer endpoint.

## Why 4000M is premature on the bound side

With `y=sqrt(Delta) R`, the exact radial equation is `y''+Q y=0`, where

\[
Q=\frac{K^2+1-a^2}{\Delta^2}
-\frac{\mu^2r^2+A+a^2\omega^2-2am\omega}{\Delta}.
\]

Its long-range form is

\[
Q=k^2+\frac{2\beta}{r}+\frac{C_2}{r^2}+O(r^{-3}),\quad
\beta=2\omega^2-\mu^2,\quad
C_2=12\omega^2-4\mu^2-a^2k^2-A.
\]

At orbit42.1 for scalar22, `k²=-3.459379979e-5`, `beta=.0899308124004`, `C2=-5.280392312`. The approximate outer turning point is `5169.718789M`; using the exact Q gives `5169.705640M`.

| Boundary radius | Exact Q | Region |
|---:|---:|---|
| 4000M | +1.004139300e-5 | Oscillatory |
| 8000M | -1.219362643e-5 | Evanescent |

Thus, even though the *asymptotic* solution decays at infinity, `4000M` lies on the oscillatory side of its outer turning point. Imposing a truncated decaying infinity expansion there can select a substantially different homogeneous combination. Moving to8000M changes this situation; proximity of the premature result to a published picture is not evidence of correctness.

The full18-mode field below changes only scalar22. The metric source and the remaining17 fields are identical. The error measure is `RMS_phi(|field|-reference_amplitude)/RMS_phi(reference_amplitude)` using the original raster's fixed angle/radius samples.

| Boundary at orbit42.1 | Error at r50 | Error at r150 |
|---|---:|---:|
| Coulomb32k | 63.14% | 59.53% |
| Leading, Rout4k | 13.76% | 13.64% |
| Li exact-prefactor order4, Rout4k | 20.67% | 14.97% |
| Li exact-prefactor order4, Rout8k | 63.14% | 59.53% |
| Li exact-prefactor order4, Rout32k | 63.14% | 59.53% |

The four-order bound field at8k differs from Coulomb32k by at most `3.5e-7` in modulus at the prescribed points; at32k by less than `4e-10`. The formal series last term alone is not a complete inner-field error estimator: in the forbidden region the selected decaying solution can become insensitive to the finite-end logarithmic derivative even when individual asymptotic terms still cancel strongly. Direct endpoint comparisons remain necessary.

## Why a nominal infinity coefficient can misrepresent the flux

On the propagating orbit41.1, `beta/(k² Rout)` is `1.9785853`, `.4946463`, and `.2473232` at Rout1k/4k/8k. These are not small expansion parameters.

For the leading boundary

\[
R=Z_I e^{ikr}r^{-1+i\beta/k}e^{-2ik\log2},
\]

its conserved radial charge current at the imposed endpoint is

\[
N=2\Delta\operatorname{Im}(R^*R')
=2k|Z_I|^2\underbrace{\frac{\Delta}{R_{\rm out}^2}
\left(1+\frac{\beta}{k^2R_{\rm out}}\right)}_{\mathcal F}.
\]

The usual unit-infinity-amplitude formula uses only `2k|Z_I|²`. It is justified only when that unit-amplitude normalization has actually been attained. At finite cutoff the leading relative correction is `(beta/k²-2)/Rout`; near threshold its coefficient diverges as `k^-2`.

At Rout4000 the exact analytic factor is `F=1.493899075774`; the integrated source-free current gives `1.493899075543`, agreement to `1.6e-10` relatively. The0th-order4000 result's apparent flux agreement is therefore partly a normalization artifact.

| Orbit41.1 boundary | Nominal total flux / Li | Conserved-current total / Li |
|---|---:|---:|
| Coulomb32k | 1.218203 | 1.218203 |
| Leading Rout4k | 1.030504 | 1.193072 |
| Li exact-prefactor order4 Rout4k | 0.747076 | 1.199477 |
| Li exact-prefactor order4 Rout8k | 1.154702 | 1.218877 |
| Li exact-prefactor order4 Rout32k | 1.218062 | 1.218204 |

Only scalar22 was replaced in these totals; the other five radiating modes remain fixed. The current is checked at two radii outside the source,400 and800M. For deliberately premature conditions, even the conserved current is the current of the chosen finite-boundary solution and is not automatically the retarded infinity answer.

At32k the paper's4th-order propagating field differs from the Coulomb field by at most `3.1e-6` in modulus at the prescribed points; the nominal total flux differs by about `.0116%`. This supports convergence toward the existing Coulomb result, not toward the reference discrepancy. No outer radius was tuned to match the reference.

## Numerical controls and limits

The16 initial cases have Wronskian relative variation at most `1.42e-10`; continuous versus original-Gauss ZI differs by at most `2.6e-11`, and ZH by `8.0e-10`. Such checks certify the integration of each chosen boundary condition, not that the condition represents infinity. Complex R and R' at20,50,100,150,200M, ZI, ZH, and both flux definitions are in the JSON. The existing calibrated reference has no r200 ring, so no reference error is invented there.

The remaining Li comparison may include source-truncation, gauge/source-convention, or background differences. These experiments do not establish which one the original authors used, and do not justify changing production to an unconverged boundary.

Artifacts: `scalar22_outer_boundary_control.json`, `scalar22_Li_order4_boundary_control.json`, `scalar22_outer_boundary_interpretation.json`, and the corresponding PNG/PDF in `docs/field_alignment_20260917/`. Generators are `src/report_li_scalar22_outer_boundary_control.py`, `src/report_li_order4_boundary_control.py`, and `src/report_li_outer_boundary_interpretation.py`.
