# Background-specific finite flux inventories

report_flux_coverage.py now accepts --alpha 0.2/0.3 and --background threshold-kerr/schwarzschild-frozen. Defaults preserve the original alpha=0.3 threshold-Kerr selection. Each inventory records a, the temporal cloud frequency and, for frozen Schwarzschild, the complex eigenfrequency, normalization and nonzero-KG-defect approximation. Source reports must match the selected background metadata and satisfy omega=omega_c+(m-1)Omega_p; incompatible channels are excluded or rejected rather than summed.

The same production discretization remains required: nr8, h32, logarithmic first panel, inner offset 0.0005M, outer source 320M, cloud mass 1, explicit metric and angular orders, explicit per-m Green settings. This change does not promote earlier L4/h16 Schwarzschild runs to production data.

Validation: the default r0=10M, L18, nt18 inventory reproduces both existing finite totals exactly: infinity 7.332063442091086e-6 and horizon -6.357677060534488e-5. For r0=20M the alpha=0.2 Kerr and alpha=0.2/0.3 frozen-Schwarzschild production inventories each correctly report 0/9 infinity and 0/18 horizon channels, with null totals. These are missing-data inventories, not newly computed physics or evidence of paper agreement.

The --output option permits regression or temporary inventories without overwriting existing snapshots. Nondefault alpha/background get distinct default filenames. Figure 2/3 still require full scans and convergence; frozen Schwarzschild retains the explicit limitation of non-exact stationary balance.
