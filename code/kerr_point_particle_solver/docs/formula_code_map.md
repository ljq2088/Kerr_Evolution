# C0 formula-to-code map

- Status: owner accepted；C0 公式核、C1/C2 实现与 plan-to-code 审计已完成
- Scope: formula/API preparation before any time loop

| Object | Findings authority | Candidate code | Contract check |
|---|---|---|---|
| Fixed defaults, internal `M=L=1`, `q_mass=1e-5`, cache key | numerical PLAN Sections 2 and 6; target Ori--Thorne generator; A6 post-ISCO interface | `config.py` | config round trip, normalization, stable key |
| `Delta`, `Sigma`, horizons, `r <-> R` | `CONVENTIONS.md`, Background and Coordinates | `geometry.py` | inverse radius check |
| Prograde ISCO quantities | A6, Post-ISCO initial-data interface | `geometry.isco_quantities` | Schwarzschild closed values; transition generator remains external |
| Inward equatorial BL and regular hyperboloidal velocities | `DERIVATION.md` Section III; Appendix E.4 | `geometry.radial_state` | BL timelike norm, independent Jacobian identities, and finite-difference radial-jet checks |
| Plunge inequalities and transition provenance | A6, Post-ISCO initial-data interface | `initial_data.py` | validation and cache protocol |
| Unique `T`-trajectory rates and all A5 stage-time jets | Appendix E.4--E.5 | `trajectory.stage_sample` | rate identities and particle-mass scaling |
| Interior first/second derivative matrices | A6 constrained two-field contract; discretization is W02.3 engineering | `operators.py` | analytic polynomial derivatives |
| Boundary-row extension point | A6 Characteristics/Axis sections | `operators.BoundaryRowExtension` | deliberately raises `OpenBoundaryError` |
| `A,B,C,D,E_R,F` | A6, Complete second-order PDE | `field.coefficients` | independent expressions at an interior point |
| Angular operator and `P/Q` identities | A6 PDE and continuous first-order system | `field.py`, `evolution.py` | independent angular expression and reconstruction of `psi_T` |
| Gaussian-dressed amplitude/phase jets | A5 Common definitions; Appendix E.7 | `source.gaussian_dressed_jets` | used by four-block and mass-linearity checks |
| LL/LJ/JJ/JL four-block source | A5 Blocks 1--4 and total source | `source.evaluate_source_blocks` | separate blocks, exact sum interface, finite interior core, unit-`mu` linearity |
| SCRI+ slice and fixed-`(ell,m)` projection | A6 SCRI+ projection | `extraction.py` | independent closed `_(-2)Y_22`, `_(-2)Y_44`, target-mode recovery, and adjacent-ell leakage |

The historical `NewPointParticlesm2Evolution` repository was not used to define the scientific formulas in this implementation.

保留的证据边界见正式 findings：A1 signature/Weyl-scalar 与 sourced-unit bridge 仍为
`conditional-open`；owner 已取消其余原定数值 validation gates，而不是将未执行检查记为通过。
