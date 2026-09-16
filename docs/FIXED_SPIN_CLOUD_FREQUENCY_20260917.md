# Fixed-spin cloud frequency and threshold radial scales

At fixed `a/M=0.88`, `M mu=0.3`, the scalar `|211>` quasibound frequency from the generalized Leaver minimal-solution condition is

`M omega_c = 0.29629353472561146 + 2.21660939645e-9 i`.

The convention is `exp(-i omega t)`: the positive imaginary part is superradiant growth. The display precision must be distinguished from finite continued-fraction truncation. The calculation used 65-digit arithmetic and a ten-component complex angular eigenproblem.

| Radial continued-fraction terms | Re(M omega_c) | Im(M omega_c) |
|---:|---:|---:|
| 150 | 0.29629353472811226847 | 2.21661280224737e-9 |
| 300 | 0.29629353472561595678 | 2.21660940302403e-9 |
| 500 | 0.29629353472561146281 | 2.21660939644579e-9 |

The last two changes are `2.49631e-12` and `4.49396e-15` in the real part, and `3.39922e-15` and `6.57823e-18` in the imaginary part. Root residuals near `1e-64` refer to the finite truncated equation, not the infinite recurrence. An independently seeded two-sided shooting check has also completed, without using this Leaver root as input. Its refined result is `0.29629353472561143 + 2.216609396778975e-9 i`, differing from the N=500 decimal output by about `3.28e-17` in the real part and `3.33e-19` in the imaginary part. The shooting refinement changed both radial endpoints, matching radius, error tolerances, boundary orders, and angular basis. This independent result supports the selected quasibound root; its evidence is `docs/environment_reproduction/fixed_spin_kerr_cloud_shooting_20260917.json`. It reuses primitive radial coefficients/boundary series but not the Leaver recurrence.

The existing `paper_leaver_cloud.py` only solves the real synchronous threshold condition, so it cannot directly accept fixed spin. The new diagnostic uses the general Frobenius reduction already implemented in `paper_leaver_radial.py`, together with the backwards minimal-solution recurrence in `paper_leaver_cloud.py`. It diagonalizes the complex angular matrix, checks that the factored radial potential is linear in the compact coordinate, and recovers the specialized real synchronous continued-fraction residual at its independently computed synchronous root. Production files were not modified.

## Comparison with the stationary synchronized cloud

The independent synchronous solution has `a/M=0.87715302759492745`, `M omega_c=0.29629324847975069`; this agrees with the production values to about `1e-14`. Raising the spin to `0.88` increases the real cloud frequency by `2.86245860768e-7`. The six-digit display `0.296294` lies another `4.65274388537e-7` above the actual fixed-spin real frequency. It must not be treated as an exact spectral input.

For scalar `m=2`, metric `m_g=1`, use

`omega = Re(omega_c) + Omega_p`, `Omega_p = 1/(r_p^(3/2)+a)`,

`kappa = sqrt(mu^2-omega^2)`, `beta = 2 omega^2-mu^2`, `eta = beta/kappa`.

Both the cloud frequency and orbital frequency change when the spin changes. The table includes both changes:

| Orbit | Synchronous kappa | Fixed-a=0.88 kappa | Relative change | Synchronous eta | Fixed-a=0.88 eta | Relative change |
|---:|---:|---:|---:|---:|---:|---:|
| 42.1 | 0.005881649411 | 0.005868971506 | -0.215550% | 15.29006680 | 15.32314653 | +0.216348% |
| 41.8 | 0.003325697165 | 0.003303295404 | -0.673596% | 27.05534359 | 27.23891311 | +0.678496% |

Using the rounded cloud frequency instead of the fixed-spin eigenfrequency produces *additional* eta errors of `+0.408260%` at `42.1` and `+1.304821%` at `41.8`. The JSON separately records a control that changes only the cloud frequency, and the complex kappa if the cloud growth rate is retained.

These numbers establish the exact spectral input and the sensitivity of asymptotic radial parameters. They do not by themselves quantify a field-amplitude or node-position correction: no forced Green solve, metric source, cloud normalization, or field reconstruction was recomputed here. Near a response resonance a small eta shift can have a larger effect, which requires an actual controlled field solve rather than an assumed scaling.

Artifacts: `src/report_fixed_spin_cloud_frequency.py` and `docs/environment_reproduction/fixed_spin_cloud_frequency_20260917.json`. The JSON preserves long decimal output, radial truncation comparisons, synchronous regression, and implementation hashes.
