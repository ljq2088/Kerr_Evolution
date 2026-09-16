# Fresh nonstatic channels for the Li spatial-field comparison

Two new orbital radii, `r_p/M=41.1` and `42.1`, have been calculated with the current production source and radial solver. Each orbit has all 16 nonstatic equatorially allowed scalar channels with `2≤ell≤5`; the static `(3,1)` and `(5,1)` channels are provided separately. No amplitude rescaling, phase fit, mode deletion, or change of the radial equation was used.

The cloud is the existing stationary synchronized alpha=0.3 solution: `a/M=0.8771530275949366`, `omega_c M=0.29629324847975713`. This differs from Li's figure label `a/M=0.88`. The calculation therefore aligns the orbital radii and stated angular sum, while retaining this explicit physical-parameter difference.

## Fixed finite resolution

Metric ellmax 6; 12-point angular projection; 8-point Gauss radial panels with a 32-point logarithmic horizon panel; source support from `r_++0.0005` to `320M`; 96 total source nodes per scalar channel. The source grid is split exactly at the particle orbit. Scalar `m=2` uses an outer boundary at `32000M` and the existing Coulomb/Whittaker boundary. All other channels use `4000M` with the inverse-radius boundary series. Every series channel passed the production last-term diagnostic; none required fallback. This is not a cutoff-convergence proof.

The two orbital processes used four metric-radius workers each. Each positive metric Fourier mode `m_g=1..6` was sampled only once, then reused for its scalar projections and its conjugate negative metric mode. Negative scalar-m responses were projected and solved separately; the complex cloud was not conjugated. Another process precomputed the final `m_g=5,6` metric caches before these were read; its manifests are retained next to the execution manifests. Execution completed successfully in about 32.8 minutes.

## Validation and scope

All 32 reports have the current recorded production source fingerprint, 96 finite complex source samples, matching saved source/file SHA-256 values, and the complete expected mode set. The largest sampled Wronskian relative variation is `1.22e-9`. The six-positive-m subset is `(2,2),(3,3),(4,2),(4,4),(5,3),(5,5)`; it retains `(2,2),(4,2)` above the continuum threshold as required for the threshold-side comparison.

| Orbit | Positive-six horizon flux | Positive-six infinity flux | All-16 nonstatic horizon flux |
|---:|---:|---:|---:|
| 41.1 | -2.46656312121e-6 | 3.26663029208e-5 | -2.71756172839e-6 |
| 42.1 | -2.50479277590e-6 | 1.83425591453e-5 | -2.75956223411e-6 |

Fluxes retain the production normalization per `q² (M_cloud/M)`. These are fixed-resolution outputs, not certified converged fluxes or a completed validation of the Li figures. The negative/zero scalar-m channels have no infinity flux here but contribute to the local field and horizon flux. Static contributions are outside this nonstatic report. The historical L6/L18 data changed several controls simultaneously and are not used as a convergence bound for these new runs.

## Artifacts

- Orchestrator: `src/report_li_field_alignment_runs.py`.
- All scalar response/source reports: `docs/field_alignment_20260917/nonstatic/`.
- Completed execution/physical-parameter/cache manifests: `nonstatic_execution_manifest.json`, `nonstatic_rp41p1_manifest.json`, `nonstatic_rp42p1_manifest.json` in the parent directory.
- Logs: `nonstatic_rp41p1.log`, `nonstatic_rp42p1.log`.
- Final independent file-integrity/coverage summary: `nonstatic_validation.json`.
- Separate metric-cache warmers' evidence: `metric_warm_rp41p1_manifest.json`, `metric_warm_rp42p1_manifest.json` (owned by the radial audit).

The diagnostic code does not modify production physics modules. Source provenance and original calculation fingerprints are preserved.
