# Lorenz reconstruction conditioning audit, 2026-09-15

This diagnostic evaluates 30 local points at rp=20, mg=1, a=0 and a=0.8771530275949366; individual ell=1,2,18; r=r_plus+0.001,3,19.9995,20.0005,40; theta=1.1. It splits the production formula into scalar CKY/kappa/compact-chi, full-current spin1, and spin2 AAB/DKW-vector/chi contributions. No production algorithm is modified.

Every covariant metric contribution is transformed to an orthonormal ZAMO frame before taking a Frobenius norm. The frame identity E g E^T=diag(-1,1,1,1) is checked. The terminal-sum condition is sum_j ||h_j|| / ||sum_j h_j||. Component ratios exclude components smaller than 1e-6 of the total norm; near-zero components cannot be used as evidence of a physical relative error.

## A real high-ell conditioning problem

At the Kerr orbit immediately inside/outside, ell18 has terminal-sum condition 2.56058e9 / 2.28420e9. Spin2 alone has condition 3.47522e6 / 3.37730e6, then the spin sectors cancel by another factor around 1.4e3. The high condition at r=3 (1.52751e10) corresponds to an extremely small ell18 metric norm, 3.25e-21, and does not imply a large physical field there.

By comparison, the dominant ell1 degree has maximum condition about 560 in these samples. Kerr ell2 reaches 5.24e4. Schwarzschild ell18 also has conditions around 2.3e9–2.5e9 near the orbit. Thus this conditioning defect is not intrinsically unique to Kerr.

## Finite sensitivity tests

The finite sensitivity tests below call the unmodified production nonstatic_metric directly. To probe internal differential operations as well as final summation, homogeneous_field_jet inputs R and R' are perturbed independently by complex relative errors of magnitude at most sqrt(2)*1e-12, using three deterministic seeds. Sourced amplitudes, trace/kappa radial data and all angular inputs remain fixed. This deliberately probes sensitivity to inconsistent radial input errors; it is not an estimate of the true error of a valid MST/GSN/TEUK solve.

| Case | Relative local metric norm changes from radial perturbation |
|---|---:|
| Kerr ell1, immediately inside orbit | 2.29e-11 to 3.34e-11 |
| Kerr ell1, immediately outside orbit | 9.58e-11 to 2.46e-10 |
| Kerr ell18, immediately inside orbit | 0.268 to 0.363 |
| Kerr ell18, immediately outside orbit | 0.246 to 0.401 |

Schwarzschild ell18 shows comparable amplification. The sensitivity can exceed the terminal-sum condition because each large contribution itself contains high-order differential operations and cancellations. A radial ODE residual of 1e-12 does not establish a reconstructed high-ell metric accuracy of 1e-12.

A second test changes Taylor coefficient arithmetic from complex128 to complex long double while retaining the same double-precision radial/angular inputs. Kerr ell18 changes by 2.64229e-3 inside and 2.11345e-3 outside; ell1 changes only 1.52e-12 / 6.65e-12. The algebraically equivalent split-term diagnostic differs from production by 4.21156e-3 / 3.50581e-3 at these points, showing that even changing summation/differentiation order matters. The sensitivity table and arithmetic comparison therefore use production calls directly, not the repartitioned expression. This is a directly observed precision dependence of the production reconstruction at high ell. It does not validate long double as the final accurate answer, because radial/angular inputs remain double and condition numbers are large.

## What it explains and what it does not

This supplies a concrete algorithmic cause for high-ell numerical contamination: large spin2 reconstruction terms and different spin sectors nearly cancel, and independently computed radial data are not protected against amplification when derivatives are composed. It supports using consistent related-spin radial data, exact Teukolsky-Starobinsky identities, analytic cancellation before floating-point evaluation, or the independently checked explicit tetrad reconstruction, with precision carried consistently from inputs through reconstruction.

It does not establish the root of the approximately 32.66% total horizon flux mismatch. The tests concern individual local ell contributions, not the final summed source integral. The dominant metric dipole has much milder conditioning, and the Schwarzschild high-ell path is similarly ill-conditioned. The actual source projection and horizon response must be recomputed under an independently stable algorithm to claim a causal correction of the main discrepancy.

Data: reconstruction_conditioning_audit.json. Reproduce with `PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_reconstruction_conditioning_audit.py`.

## Actual dipole radial-backend propagation cross-check

The companion actual 88-node source integration `metric_backend_response_L1_s2AUTO_g*.json` changes the gauge-sector radial backend while retaining the L1 physical problem. The AUTO result is FH=-6.805625878766075e-5, ZH=-0.05003928029885452+0.03370367360681479i; direct TEUK gives FH=-6.805625878772487e-5, ZH=-0.05003928029887946+0.03370367360682864i. The horizon flux change is about 9.4e-13. This independently argues against HBL versus direct-TEUK radial integration as the main dipole discrepancy. Shared boundary normalization and theoretical source construction remain outside that comparison. Consequently the high-ell sensitivity documented here is a subsequent precision-repair item, not an established cause of the 32% total discrepancy.
