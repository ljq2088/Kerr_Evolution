# Three-orbit finite Figure 2 comparison

report_three_orbit_flux.py independently resums the complete source-file fluxes selected by the current 9-channel infinity / 18-channel horizon coverage reports at r0=10,20,30M. It verifies the coverage sums against the actual completed files, and records all source fingerprints. The common paper reference is paper_figure2_selected_radii.json; no amplitude fit is applied.

| r0/M | Infinity magnitude excess | Horizon magnitude excess |
|---:|---:|---:|
| 10 | 4.5691% | 40.6688% |
| 20 | 5.0454% | 32.6552% |
| 30 | 5.2437% | 32.0855% |

The reference gives horizon magnitude only; computed signed values remain in figure2_three_orbit_comparison.json. This version consistently uses direct Figure 2 readouts at all three radii. Small differences from older quoted r20 percentages reflect previously used reference readouts, not a rescaling of the computed flux.

figure2_three_orbit_comparison.png was rendered and visually checked. It shows only three sampled orbits, not a complete radial scan. The digitized paper values have no rigorous readout error bound; source/metric convergence and the discrepancies remain unresolved.

The earlier rp30_infinity_reference_comparison.json is a historical snapshot captured before its final horizon mode completed. Its coverage SHA256 has been checked against git revision 740afd3 and now explicitly points to that revision. The new comparison uses the current complete coverage.
