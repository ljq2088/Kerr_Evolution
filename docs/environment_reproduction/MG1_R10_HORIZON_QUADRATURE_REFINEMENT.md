# r0=10M near-horizon quadrature refinement

This is a completed six-channel numerical comparison, not global convergence or agreement with the paper. The logarithmic horizon panel changes from 32 to 64 quadrature nodes (96 to 128 total source nodes). Other recorded physical and numerical parameters agree. However, the new run also uses the shared angular-basis extension described in ANGULAR_BASIS_EXTENSION.md, so this is not a strictly isolated quadrature experiment.

| (ell,m) | Infinity relative change | Horizon relative change | Common normal-panel source relative L2 |
|---|---:|---:|---:|
| (2,2) | -2.107133074e-07 | 1.966420309e-07 | 2.708664505e-07 |
| (4,2) | -2.173180458e-08 | 3.667690824e-06 | 8.765228124e-07 |
| (6,2) | 6.293504162e-09 | 1.790561854e-05 | 1.465339948e-06 |
| (0,0) | zero in both runs | 1.442175046e-08 | 4.626394849e-07 |
| (2,0) | zero in both runs | 8.601974688e-07 | 7.072595828e-07 |
| (4,0) | zero in both runs | 7.423302518e-06 | 4.915354906e-07 |

The dominant (0,0) horizon contribution changes from -6.272275267097329e-5 to -6.272275357554517e-5, a relative change of 1.442175046e-8. This change does not explain the approximately 40.7% excess of the existing finite horizon sum over the digitized paper reference at this orbit. Weak-channel relative changes should be interpreted with their absolute fluxes, retained in the JSON report.

The 64 unchanged normal-panel radii were checked explicitly; their recomputed sources are close but not identical. Thus the tiny flux difference must not be interpreted as a rigorous isolated quadrature error bound. The horizon inner cutoff, source convention, metric reconstruction, normalization, and complete mode sum still require independent checks.

Inputs and SHA256 fingerprints are in mg1_r10_horizon_quadrature_refinement.json. The completed batch manifest is scalar_batch_mg1_604b085df0b7.json. No fit or normalization adjustment was applied.
