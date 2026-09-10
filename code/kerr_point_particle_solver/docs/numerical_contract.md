# Numerical contract

- Status: owner accepted；C0--C2 已实现，原计划中未执行的其余数值 gates 由 owner 取消
- Runtime arrays: NumPy `float64` coordinates and `complex128` fields/sources
- Shape: field arrays use `(n_R, n_y)`; radial derivatives use axis 0 and angular derivatives axis 1
- Units: runtime geometry is dimensionless with internal `M=1`, `L=1`; `mass_scale` is retained only for later physical-scale restoration
- Fourier convention: `psi_4 = sum_m psi_{4,m} exp(i m Phi)` and source phase `exp(-i m Phi_p)`
- Trajectory: only `(R_p(T), Phi_p(T))` is dynamical; `r_p=L^2/R_p` is a derived evaluation coordinate
- Source amplitude: linear in unit particle mass `mu`; it is separate from the fixed transition mass ratio
- Transition prescription: fixed to `ori-thorne-v1` with `q_mass=1e-5`; neither is a run-time input. A versioned flux model is required before an initial-data cache key or manifest is accepted
- Post-light-ring duration: default fixed to `120 M` by the owner and implemented in `EvolutionConfig`; callers may record an explicit override only in later reviewed numerical work
- Endpoint: `T_bound=T_LR^scri+DeltaTpost`, `N_end=floor(T_bound/DeltaT_out)`, and `T1=N_end DeltaT_out<=T_bound`. The default driver stops at `T1`; `T_bound,T1,N_end,DeltaTpost` are recorded in run and mode metadata.

## C0 公式核

C0 固定了 Kerr 几何、轨迹变量、A--F 系数、四块 source、$P/Q$ 一阶化和
$\mathscr I^+$ 投影的公式到代码映射。其独立审阅见
[`formula_code_map.md`](formula_code_map.md)。C0 曾保留的 boundary-row extension point
已在 C1 中由实际 production rows 替代。

## C1 数值实现

- Target-consistent Ori--Thorne inputs are fixed to `q_mass=1e-5`, `lmax=6`, `T_OT=3.412`, and `X_start=-0.5381356808642278`. Cache keys include these values and all flux/schema/generator versions.
- Trajectory uses SciPy DOP853 with `rtol=2e-11`, `atol=2e-13`, `max_step=0.25M`; events are light ring, horizon/source-off, and `R_H+sqrt(2 log(1e14))*sigma_R` terminal.
- Spatial grids use seven points. Ordinary weights use the locally scaled Fornberg recursion; factor-aware rows solve a locally scaled seven-function basis and record its condition number. All matrices are sparse CSR with at most seven entries per row. Radial rows are centered/shifted polynomial rows including `R_H` and `R=0`; the three axis-adjacent field rows on each side use the corresponding local spin-weight factor; source angular rows remain ordinary shifted rows. D1 and D2 are independent.
- Source startup defaults to `tau_on=20M` with hard API range `0<tau_on<=40M`; the same C^2 quintic is applied at amplitude-jet level. Gaussian support uses the two-dimensional `epsilon_G=1e-14` ellipse. SCRI source is analytic zero when retained support does not reach it. The dedicated horizon-total evaluator directly transcribes all four regular A5 operands and their sum at `R=R_H`, with `rho^8/R^4` evaluated as `R^4/(L^2+i a R y)^8`; it does not consume generic horizon-block output. Particle-center `R_p>=R_H` returns four exact zero blocks.
- Classical RK4 resamples the trajectory/source and recomputes spatial derivatives at every stage. The global timestep evaluates motion/phase/base-amplitude extrema only on the fully-on interval `[T_stable,T_H)`. `T_stable=100M` by default; for `tau_on<T_H<120M`, it is `max(tau_on,T_H-20M)`; `T_H<=tau_on` is rejected. Extrema use all DOP853 segments, endpoint/midpoint probes, bounded optimization on the 16 largest candidates, and a 1.01 envelope. The fifth limit is `dt_ramp=tau_on/40`; startup-window derivatives do not enter global `Gamma_S`. Other factors remain `C_CFL=0.20`, `f_R=0.50`, `f_Phi=0.10`, `f_S=0.15`, with `Q_floor=1e-12 Q_ref`.
- Science mode sampling defaults to `0.1M`, full-field memmap cadence to `5M`, and checkpoint cadence to `10M`. Axis endpoints are reconstructed before augmented-uniform-grid Simpson projection. Disk writes occur only through accepted-level callers, never inside `rk4_step`.
- Current target modes are restricted to `m=2,4`. Science mode/SCRI samples remain on the exact uniform output grid; a non-grid final state is stored only in run metadata and the final checkpoint. Checkpoints persist separate mode, SCRI, and full-field valid counts for restart truncation.
- Mode, SCRI, run, and checkpoint schemas carry completion state and provenance hashes. Rolling checkpoints retain two atomic generations. Restart truncates accepted outputs to the checkpoint index and preserves mode/SCRI/full-field history.
- The code hash is SHA-256 over the actual sorted package source filenames and contents; restart verifies it. It is not a fixed version label.
- Run metadata records the actual Python, NumPy, and SciPy versions. Read-only mode, SCRI, and full-field loaders validate schemas, completion state, hashes, dtype/shape, and full-field valid-slice counts; non-complete test data require an explicit loader opt-in.
- Artifact status is the closed enum `incomplete|test|complete`; unknown values are always rejected. Restart also requires exact Python/NumPy/SciPy environment equality. Mode, SCRI, and full-field files shorter than checkpoint counts are rejected; longer files are atomically truncated and revalidated to the exact checkpoint counts before evolution resumes.
- Every persisted count is validated before integer conversion: it must be a scalar real value, non-boolean, finite, exactly integer-valued, and nonnegative. Mode/SCRI files persist `sample_count`; checkpoints persist output/mode/SCRI/full-field counts; full-field sidecars persist `valid_slices`; run metadata persists endpoint/output counts. Fractional, negative, string, and boolean values are rejected.

## C2 strain 实现

- Output is asymptotic `H_lm=r h_lm`, not finite-distance `h`; the explicit bridge is `H_ddot=L^2 psi4`.
- FFI uses the mode file's actual `T_a=T_stable` (default `100M`, approved fallback as above), `omega_i=|m Phi_dot(T_a)|`, baseline `omega0=0.75 omega_i`, a two-sided `10M` quintic taper, and conditional-open A1 metadata plus an input hash. Conversion rejects data that do not cover the recorded analysis start.
- C2 sets `T_b=T1` and records analysis bounds, FFT convention, taper/cutoff, sampling status/results, input schema/hash, source events, normalization, and source error metadata. A short/test mode is rejected by default; explicit test conversion remains marked `test`. Only a mode from a complete full run may produce `status=complete`.
- The sampling comparison separately integrates `0.05M` data and its every-other-point `0.10M` restriction. Its threshold is `min(1e-3,0.1 epsilon_V2)` and therefore requires a positive V2 estimate.
- The coarse sampling grid is selected by absolute times `T0+n*0.10M` with absolute tolerance `2e-11M`; it is not `T_fine[::2]`. Fine spacing must be uniformly `0.05M`, and the selected coarse times must be uniformly `0.10M`. Thus a driver series beginning at `0.05M` selects `0.10,0.20,...`.

## 证据边界

SXS:BBH:0305 的精确 remnant spin 已执行 $\ell_{\max}=6$ BHPT ISCO flux，
$(2,2)$ 与 $(4,4)$ medium production runs 均已完成。Owner 以最低分辨率跨代码
波形一致性接受 solver，并取消原计划中其余 V1--V3 数值 gates；这些未执行检查不
记为通过。A1 strain sign/normalization bridge 仍在 C2 metadata 中保持
`conditional-open`。
