# Massive scalar radial precision audit, 2026-09-15

This audit holds the 21 stored source arrays and finite quadrature fixed at alpha=0.3, a=0.8771530275949366, rp=20. It changes the homogeneous massive scalar solver only. The source files are hashed and are not relabeled as freshly reconstructed sources. The earlier independent Leaver proof for the dominant scalar (0,0) is not repeated.

## Equation and algorithm regime

The metric reconstruction uses massless spin-weighted Teukolsky modes at omega_g=m_g Omega_p. The forced environment has a different frequency omega_s=omega_cloud+m_g Omega_p, scalar mass mu=0.3, and asymptotic wave number k=sqrt(omega_s^2-mu^2). Massless MST/GSN cannot be substituted unchanged into this massive Klein-Gordon radial equation.

The current environment solver uses DOP853, a fourth-order ingoing Frobenius expansion, and a sixth-order inverse-radius infinity expansion. A Coulomb/Whittaker alternative is already available near the mass threshold. At the tested point, min |k|=0.0668823504 and min |omega_s-m_s Omega_H|=0.2852214879 across all 21 channels. Neither the mass threshold nor horizon synchronization is approached by these nonstatic channels.

For contrast, the scalar (2,2) radiation threshold is rp=41.6608484. At rp=41.6 the far-expansion parameter |(2 omega_s^2-mu^2)/(k^2 rmax)| at rmax=1000 is 18.50 and the retained last-term ratio is 1.004; at rp=41.8 these are 8.14 and 1.089. The existing producer's unresolved-series guard rejects such runs. At rp=20 these are 0.02212 and 5.44e-10 for this channel. These threshold failures cannot explain the rp=20 discrepancy.

## Fixed-source tests

Relative changes below are dimensionless ratios, not percent. Totals use the same boundary-dependent mode coverage as the original paper comparison (horizon ell<=5, infinity ell<=6).

| Modification | Relative horizon total change | Relative infinity total change |
|---|---:|---:|
| rtol 1e-11 to 5e-14 | -4.82309e-11 | -5.10738e-11 |
| rmax 1000 to 2000; inner offset 1e-4 to 1e-5; infinity order 6 to 12; rtol=1e-12 | -7.97167e-11 | -3.53026e-11 |
| Nine propagating channels use Coulomb/Whittaker at rmax=4000 | +1.50768e-12 | -3.58914e-9 |

Baseline horizon total=-6.866732369763171e-5; infinity total=1.9210336966212187e-5. The last row retains baseline bound-channel contributions. The Coulomb approximation keeps the far potential only through r^-2, so it is a cross-check, not an exact replacement.

For every individual channel, tightening rtol changes the horizon flux by at most 1.68e-10 and infinity by 1.14e-10. Boundary refinement changes horizon by at most 1.37e-10 and infinity by 1.34e-9. No substantial algorithm-dependent change appears.

The maximum Wronskian subtraction condition (|u v'|+|v u'|)/|u v'-v u'| over actual source nodes is 2.172. The dominant scalar (0,0) horizon integral has cancellation condition sum |w_i J_i Up_i/W| / |Z_H|=1.617. Across all modes it is at most 3.851. Large absolute normalization of bound Up therefore does not correspond to a catastrophic cancellation at this test point. Baseline Wronskian spread is at most 1.31e-9 and drops below 1e-11 when tolerance is tightened.

## Limits

These tests exclude the specific scalar radial precision and endpoint choices as an explanation of the 32.66% horizon disagreement at this point. They do not independently validate the metric, source samples, quadrature, finite mode cutoff, or author conventions. Narrow resonances and k=0 require additional treatment. Increasing the bound-state outer radius indefinitely is unsafe because unscaled exponentially large homogeneous solutions can overflow; this audit uses rmax=2000 for bound channels.

The JSON separately records a sum of all 21 channels at both boundaries. Its horizon total includes the additional weak scalar (6,2) contribution -1.66673e-11, absent from the original horizon ell<=5 coverage. This coverage difference is explicit and does not materially alter the main discrepancy.

Reproduce with `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_massive_radial_precision_audit.py`.

Data: `massive_radial_precision_audit.json`. No production module is changed by this diagnostic.
