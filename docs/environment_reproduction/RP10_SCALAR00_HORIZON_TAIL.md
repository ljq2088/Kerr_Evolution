# r0=10M scalar (0,0) omitted horizon integral audit

The existing two-phase Frobenius extrapolation was applied to the completed L18, nt18, nr8 source data with 32 and 64 logarithmic horizon nodes. The finite inner cutoff is 0.0005M. No source or production flux file was overwritten. Input SHA256 fingerprints and all polynomial orders/windows are in rp10_scalar00_horizon_tail_audit.json.

| Horizon nodes | First-order window | Corrected horizon flux | Relative change from finite cutoff | Held-out integrand absolute residual |
|---:|---:|---:|---:|---:|
| 32 | 5 cutoff lengths | -6.275060322927e-05 | 4.440264036e-04 | 4.880174083e-06 |
| 32 | 10 cutoff lengths | -6.275064188801e-05 | 4.446427468e-04 | 1.076606870e-05 |
| 32 | 20 cutoff lengths | -6.275065752570e-05 | 4.448920614e-04 | 5.141256943e-06 |
| 64 | 5 cutoff lengths | -6.275068357721e-05 | 4.452929770e-04 | 4.958670009e-06 |
| 64 | 10 cutoff lengths | -6.275067924550e-05 | 4.452239157e-04 | 4.870153714e-06 |
| 64 | 20 cutoff lengths | -6.275067767264e-05 | 4.451988395e-04 | 4.812483395e-06 |

For the 64-node data the three first-order estimates increase the flux magnitude by approximately 0.0445%. They cluster near -6.275068e-5, compared with the uncorrected -6.27227535755e-5. This effect is much smaller than the existing approximately 40.7% finite total horizon discrepancy against the digitized paper value, and moves its magnitude upward.

The fit residuals are absolute integrand residuals, not relative flux errors. Typical integrand magnitudes are about 0.03. All orders zero, one and two and windows 5, 10 and 20 are retained, including poorly conditioned or insufficient held-out combinations. This is not a rigorous extrapolation bound and does not replace a smaller-cutoff run. The angular-basis extension also differs between the archived 32-node and new 64-node calculations, as documented in MG1_R10_HORIZON_QUADRATURE_REFINEMENT.md.

Reproduction: run src/report_horizon_tail.py on the two input paths listed in the JSON, with --output selecting a separate destination. The report tool now supports this option to preserve earlier horizon_tail_audit.json results.


A direct cutoff-halving batch has now been launched: r0=10, L18, nt18, nr8, h64, inner source offset 0.00025M, outer source 320M and Green boundaries unchanged (offset 0.0001M, outer 1000M, series). It includes scalar (2,2),(4,2),(6,2),(0,0),(2,0),(4,0), with three metric-sampling workers. The comparison baseline is scalar_batch_mg1_604b085df0b7.json at offset 0.0005M, using the same current angular-basis implementation. Preflight verified the 64 normal-panel radii and weights are exactly unchanged and the new cutoff stays outside the Green inner boundary. The 64 logarithmic-panel radii change. Log: outputs/rp10_mg1_inner0p00025_h64.log. This run is in progress; no smaller-cutoff convergence result is claimed yet.
