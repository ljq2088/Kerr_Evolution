# 模块与原作者参考资料映射（2026-09-15）

原Sam完整辐射度规代码已在实际参数a=0.8771530275949366、rp=20、m=1、L=10完成两次运行及边界精化。三处半径r=2,10,30、受截断保护的球谐j=1,2,3、十个加权投影均完成复数比较。每个投影以j=1..3组成的向量计算L2相对误差，最大为8.191e-7，不能表述为每个独立非零系数都达到相同相对误差或全度规精度1e-14。

正式粗/细边界为inford6/horord5/rinf10000与inford7/horord6/rinf20000，固定xhor=1e-4、WP50、AccuracyGoal24。精化使原作者度规变化至多3.337e-12；本地对精化结果的最大误差仍为8.191e-7。两组耗时过长的inford10符号展开实验已放弃作为完整运行候选，其此前完成的匹配输出保留有效；它们不是待办。

两套度规输入同一个本地云场和协变源公式后，截断scalar00源在r2/10/30的差异为2.478e-9、1.325e-9、1.834e-7。此处使用同一ThresholdCloud和Hessian，验证的是原作者度规输入误差的有界传播；它不是原作者环境源代码验收，且j<=3的J不是fullJ。

当前结果不支持这些已采样低阶模块因普通数值误差产生百分比差异，但并未验证完整径向网格、高j、所有m、原作者环境源或最终通量。原作者代码运行、官方依赖运行、分发样例、独立重实现与PDF数字化分别列明；当前版本也不等于原论文历史依赖锁定。

## 原代码执行状态

- `original_operator_and_matching_helpers`：completed。
- `sam_matching_a06_rp8_m1_L10`：completed。
- `sam_full_metric`：completed: actual-parameter full metric i6/h5 and refinedi7/h6 at r2/10/30,m1,L10,protectedell1..3。
- `sam_actual_a0877_rp20`：matching and both full-metric sampled comparisons completed。
- `environmental_source_and_flux_original_run`：not available。
- `official_Mathematica_BHPT_11cases`：completed: orbit, SWSH, massless MST In/Up R/Rprime and sampled Green kernel。
- `bounded_original_metric_input_source_propagation`：completed for both i6/h5 and i7/h6: same local cloud/Hessian, j1..3/scalar00; not original environment-source execution。
- `author_boundary_refinement`：completed: maximum author change3.34e-12; local-vs-fine maximum8.191e-7。
- `inford10_symbolic_expansion_experiments`：abandoned as full-metric run candidates; completed original matching outputs remain valid and recorded; not pending work。

## 已完成运行的可复核证据

- **original_scalar_kappa_helpers**（`original_author_run`）：Mathematica Kernel14 executed verbatim Sam and Conor scalar/kappa operators and ConstructSolution-style matching helpers. Scalar/kappa operator differences exactly0; three complex basis scalings through1e±25 preserve exact jumps; local matrix-form difference2.66e-16. This is operator/helper algebra, not full physical kappa integration by this helper test alone.
  结果：`docs/lorenz_public_code_manifest_20260915.json`；原始输出：`outputs/paper_metric_reference/radial_author_modules.json`；driver：`outputs/paper_metric_reference/radial_author_modules.wls`；来源：`outputs/paper_metric_reference/code_audit.json`。
- **sam_matching_a06_rp8**（`original_author_run`）：Original Sam matching phase executed at a=.6,rp8,m1,L10,WP50. Four complex spin1/kappa jumps per ell1/2/3 agree with production at row-scaled1.98e-14,1.23e-12,1.15e-13. Modes ell4..10 are excluded as cutoff-edge diagnostics; not a full metric or flux comparison.
  结果：`docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json`；原始输出：`outputs/paper_metric_reference/sam_L10_m1/matching.json`；driver：`outputs/paper_metric_reference/sam_L10_m1/run_sam.wls`；来源：`docs/environment_reproduction/sam_minimal_driver_manifest_20260915.json`。
- **official_BHPT_dependency_11cases**（`official_dependency_run`）：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.
  结果：`docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`；原始输出：`outputs/paper_metric_reference/wolfram_dependencies/official_results`；driver：`docs/environment_reproduction/wolfram/bhpt_official_probe.wls`；来源：`docs/environment_reproduction/wolfram/official_bhpt_manifest.json`。
- **sam_matching_actual_a0877_rp20**（`original_author_run`）：Original Sam matching executed at actual discrepancy a=.8771530275949366,rp20,m1,L10. Four spin1/kappa jump entries for ell1/2/3 agree with production at row-scaled errors 2.2993e-14,3.3554e-12,2.1032e-13. Six-mode cutoff guard; symmetry-zero entries use row scaling. Full metric completion is reported separately; this matching report alone has no source/flux conclusion.
  结果：`docs/environment_reproduction/sam_actual_vs_production_jumps_ell1to3.json`；原始输出：`outputs/paper_metric_reference/sam_actual_L10_m1/matching.json`；driver：`outputs/paper_metric_reference/sam_actual_L10_m1/run_sam.wls`；来源：`docs/environment_reproduction/sam_actual_minimal_driver_manifest_20260915.json`。
- **sam_actual_full_metric_i6_h5**（`original_author_run`）：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.
  结果：`docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`；原始输出：`outputs/paper_metric_reference/sam_actual_L10_i6_h5/metric_samples.json`；driver：`outputs/paper_metric_reference/sam_actual_L10_i6_h5/run_sam.wls`；来源：`outputs/paper_metric_reference/sam_actual_L10_i6_h5/driver_manifest.json`。
- **sam_actual_full_metric_i7_h6**（`original_author_run`）：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.
  结果：`docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`；原始输出：`outputs/paper_metric_reference/sam_actual_L10_i7_h6/metric_samples.json`；driver：`outputs/paper_metric_reference/sam_actual_L10_i7_h6/run_sam.wls`；来源：`outputs/paper_metric_reference/sam_actual_L10_i7_h6/driver_manifest.json`。

## 原作者数值样例与表格

- Conor固定commit的a=.6、rp8、L4样例HDF5：trace显式spheroidal→spherical后六点误差2.09e-12至9.75e-11。其他分量存在样例连续性/生成副本问题，不能据此判定论文或本地实现错误；新运行Sam代码提供单独、可重现的更强证据。
- 2024 arXiv包中的amplitudes.dat是原作者数值表；Operators.nb是原作者符号代码。已完成的算子运行使用Sam/Conor原片段，不等于执行整个Operators.nb。
- 原作者新输出的非标准JSON小数指数只在字符串外补零，例如0.e-49→0.0e-49；原始字节与SHA256保留。比较报告逐项记录语法正规化，不改数值、尺度或相位。
- 固定版本Conor仓库全部2069可读文本中未定位完整有质量环境标量流水线。旧notebook确有boson标签的度规输出路径，不能将其抹去，也不能把它们升格为完整环境程序；检索边界见专门报告。

## 逐模块清单

### Orbit and frequency bookkeeping

实现：`src/lorenz_weyl.py`, `src/environment_lorenz_mode.py`, `src/environment_source.py`。

原文：Dyson numerical-procedure appendix: omega_s=omega_c+m_g Omega_p; Kerr circular orbit in DDKW source equations。

对照输入：M,a,rp,m_g,omega_c。输出：Omega_p,u^t,E,L,omega_g,omega_s。

现有证据：
- `original_author_code` — `outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl`：OrbitalData / KerrAzimuthalFrequency; source uses M=1; no separate author orbit-array comparison yet。
- `official_dependency_run` — `docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.。

仍缺少：Exact orbital and cloud parameters used for every plotted 2025 marker.

### Cloud eigenfrequency and normalization

实现：`src/environment_cloud.py`, `src/environment_source.py`, `src/environment_schwarzschild_cloud.py`, `src/paper_leaver_cloud.py`。

原文：Dyson stationary 211 cloud and numerical appendix; epsilon=alpha^3 sqrt(Mc/M); M_cloud=-integral T^t_t Sigma dr dOmega per referenced thesis convention。

对照输入：alpha,a,n,l_c,m_c and cloud mass convention。输出：complex omega_c,radial profile,R derivatives,normalization and mass integral。

现有证据：
- `independent_reimplementation` — `docs/environment_reproduction/leaver_cloud_full_audit.json`：Local independent Leaver implementation, not author code。
- `internal_identity` — `docs/environment_reproduction/cloud_mass_measures.json`：Comparison of mass measures; not author cloud array。
- `original_author_code` — `docs/CONOR_ENVIRONMENT_SOURCE_SEARCH_20260915.md`：Complete fixed-commit readable-text search:2111files,2069texts,47.35MB,Git blob checked; no identifiable full massive cloud/Klein-Gordon environment source/flux pipeline found. Old notebooks contain real boson-tagged metric output paths, so this is not a claim of no environment-related artifacts. CompressedData/binary contents and unpublished local referenced files excluded.。

仍缺少：Raw author cloud eigenfrequency, radial profile and normalization arrays at alpha=.3,a=.877153... .

### Spin-weighted spheroidal harmonics

实现：`src/environment_angular_diagnostic.py`, `src/environment_source.py`, `src/paper_precise_angular.py`。

原文：DDKW Teukolsky separation and spin harmonics; Dyson angular spheroidal projection。

对照输入：s,ell,m,c=a*omega and harmonic phase/normalization。输出：lambda,S(theta),dS,cos/cos^2 mixing matrices and spherical-expansion coefficients。

现有证据：
- `third_party_library_run` — `docs/environment_reproduction/swsh_independent_qnm_audit.json`：Independent qnm angular implementation; not Dyson source。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_metric_aligned_comparison.json`：Explicit scalar spheroidal-to-spherical B transform makes author trace agree 2.1e-12 to9.75e-11 at six points; validates combined scalar trace/basis, not all angular sectors independently。
- `official_dependency_run` — `docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.。

仍缺少：Original production angular cutoff and complete high-ell/grid outputs. Direct official original SWSH values, slopes and eigenvalues for11bounded cases now pass, but current official package revisions are not the original paper dependency lockfile.

### Massless spin +/-2 radial and Weyl amplitudes

实现：`src/lorenz_weyl.py`, `src/lorenz_metric.py`, `src/lorenz_mode_jet.py`。

原文：DDKW eq:teukolsky-eqns,eq:TS-identities/eq:TSs; WDK2024 eq:Psi-FD and eq:TeukolskyInhomogeneousModes。

对照输入：a,rp,l,m,s,frequency,In/Up normalization。输出：R,Rprime,Wronskians,Weyl/source amplitudes,TS constants。

现有证据：
- `published_author_numeric_table` — `docs/environment_reproduction/lorenz_weyl_validation.json`：Compares original WDK2024 amplitudes.dat via documented i/omega convention bridge。
- `third_party_library_run` — `docs/environment_reproduction/radial_backend_invariants.json`：BHPT backends; not independent author environmental code。
- `official_dependency_run` — `docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_boundary_convergence.json`：Both actual original-author metric runs completed. Boundary controls6/5/rinf10000 ->7/6/rinf20000 at fixed xhor1e-4,WP50,goal24. Maximum author relative change at r2/r10/r30=6.83e-14/2.06e-14/3.34e-12, while local-vs-fine maximum remains8.191e-7. Hence the residual sampled local difference is not explained by these author boundary settings. This is not horizon-limit/high-L convergence or a final flux test.。

仍缺少：Author raw same-mode radial values and amplitudes at actual discrepancy parameters; exact MST/GSN/backend choice in 2025 adapted package.

### Massless spin +/-1 radial and vector amplitudes

实现：`src/lorenz_spin1.py`, `src/lorenz_spin1_chiral.py`, `src/lorenz_metric.py`, `src/lorenz_mode_jet.py`。

原文：DDKW eq:mp-s1,eq:TSs and coupled source jumps; WDK2024 eq:teuk-s1。

对照输入：a,rp,l,m,s,In/Up side。输出：radial values,vector source amplitudes and full chiral current。

现有证据：
- `published_author_numeric_table` — `docs/environment_reproduction/lorenz_spin1_validation.json`：Original WDK2024 table after sqrt(2)/omega^2 convention bridge。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_metric_aligned_comparison.json`：m1 HDF5 likely omits ell1 vector in a generation path; omit-only diagnostic improves some coefficients but is not a justified production change。
- `original_author_run` — `docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json`：Original Sam matching phase executed at a=.6,rp8,m1,L10,WP50. Four complex spin1/kappa jumps per ell1/2/3 agree with production at row-scaled1.98e-14,1.23e-12,1.15e-13. Modes ell4..10 are excluded as cutoff-edge diagnostics; not a full metric or flux comparison.。
- `official_dependency_run` — `docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_vs_production_jumps_ell1to3.json`：Original Sam matching executed at actual discrepancy a=.8771530275949366,rp20,m1,L10. Four spin1/kappa jump entries for ell1/2/3 agree with production at row-scaled errors 2.2993e-14,3.3554e-12,2.1032e-13. Six-mode cutoff guard; symmetry-zero entries use row scaling. Full metric completion is reported separately; this matching report alone has no source/flux conclusion.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_boundary_convergence.json`：Both actual original-author metric runs completed. Boundary controls6/5/rinf10000 ->7/6/rinf20000 at fixed xhor1e-4,WP50,goal24. Maximum author relative change at r2/r10/r30=6.83e-14/2.06e-14/3.34e-12, while local-vs-fine maximum remains8.191e-7. Hence the residual sampled local difference is not explained by these author boundary settings. This is not horizon-limit/high-L convergence or a final flux test.。

仍缺少：Full radial-domain and high-ell original vector/metric checks plus additional m. Official actual spin±1 radial samples and actual low-ell Sam jumps pass; original full spin1 metric sector at three actual radii now agrees to2e-12. Final covariant source still needs direct original data.

### Massless spin0 trace and chi

实现：`src/lorenz_trace.py`, `src/lorenz_chi.py`, `src/lorenz_mode_jet.py`。

原文：DDKW eq:radial-h,eq:hsol,eq:jumps-h; WDK2024 eq:chi-DKW and eq:teuk-s0。

对照输入：a,rp,l,m,In/Up side,frequency,angular convention。输出：h,chi and radial derivatives/source amplitudes。

现有证据：
- `published_author_numeric_table` — `docs/environment_reproduction/lorenz_trace_validation.json`：Original WDK2024 amplitudes.dat。
- `published_author_numeric_table` — `docs/environment_reproduction/lorenz_chi_validation.json`：Original table with -i/(2omega) bridge。
- `independent_reimplementation` — `docs/environment_reproduction/trace_kappa_radial_precision_audit.json`：Fixed source trace/kappa endpoint/tolerance/variation alternatives。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_metric_aligned_comparison.json`：Original 2026 HDF5 trace after documented basis conversion matches local complex trace at six points to <=9.75e-11。
- `original_author_run` — `docs/lorenz_public_code_manifest_20260915.json`：Mathematica Kernel14 executed verbatim Sam and Conor scalar/kappa operators and ConstructSolution-style matching helpers. Scalar/kappa operator differences exactly0; three complex basis scalings through1e±25 preserve exact jumps; local matrix-form difference2.66e-16. This is operator/helper algebra, not full physical kappa integration by this helper test alone.。
- `official_dependency_run` — `docs/environment_reproduction/wolfram/official_bhpt_AUTO_comparison.json`：Executed official Mathematica BHPT code, not the complete 2025 environmental code release. 11cases: a=.6,rp8,ell2,m1,s=-2..2; actual a=.8771530275949366,rp20,ell1,s0±1,m±1. Three radii/case and In/Up; 66samples. No fitted scale or phase. R max relative3.0065e-14,Rprime3.3356e-14; angular lambda absolute1.3172e-13,S5.8009e-15,dS1.3989e-14; massless Green sample relative2.1763e-14.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_boundary_convergence.json`：Both actual original-author metric runs completed. Boundary controls6/5/rinf10000 ->7/6/rinf20000 at fixed xhor1e-4,WP50,goal24. Maximum author relative change at r2/r10/r30=6.83e-14/2.06e-14/3.34e-12, while local-vs-fine maximum remains8.191e-7. Hence the residual sampled local difference is not explained by these author boundary settings. This is not horizon-limit/high-L convergence or a final flux test.。

仍缺少：Original full chi and kappa radial profiles at production parameters. Official original massless trace radial samples and author HDF5 trace now pass; massive scalar Green and full metric/source are separate checks.

### Kappa inhomogeneous resolvent and Taylor series

实现：`src/lorenz_kappa.py`, `src/environment_trace_variation.py`, `src/environment_radial_variation.py`, `src/paper_kappa_radial.py`。

原文：DDKW eq:Box-kappa,eq:kappa-lj,eq:kappa-ll,eq:radial-kappa and In/Up asymptotics; WDK2024 eq:kappa。

对照输入：a,rp,l,m,trace mode; boundary orders/radii; inhomogeneous normalization。输出：kappa and derivatives,homogeneous constants and source response。

现有证据：
- `independent_reimplementation` — `docs/environment_reproduction/paper_kappa_radial_audit.json`：Our implementation of alternative published formula, not original code。
- `internal_identity` — `docs/environment_reproduction/dipole_kappa_flux_boundary_audit.json`：Actual source propagation of boundary/analytic variation。
- `independent_reimplementation` — `docs/environment_reproduction/trace_kappa_radial_precision_audit.json`：Nine-node complex source comparison; endpoint extremely close to horizon can fail integration。
- `original_author_run` — `docs/lorenz_public_code_manifest_20260915.json`：Mathematica Kernel14 executed verbatim Sam and Conor scalar/kappa operators and ConstructSolution-style matching helpers. Scalar/kappa operator differences exactly0; three complex basis scalings through1e±25 preserve exact jumps; local matrix-form difference2.66e-16. This is operator/helper algebra, not full physical kappa integration by this helper test alone.。
- `original_author_run` — `docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json`：Original Sam matching phase executed at a=.6,rp8,m1,L10,WP50. Four complex spin1/kappa jumps per ell1/2/3 agree with production at row-scaled1.98e-14,1.23e-12,1.15e-13. Modes ell4..10 are excluded as cutoff-edge diagnostics; not a full metric or flux comparison.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_vs_production_jumps_ell1to3.json`：Original Sam matching executed at actual discrepancy a=.8771530275949366,rp20,m1,L10. Four spin1/kappa jump entries for ell1/2/3 agree with production at row-scaled errors 2.2993e-14,3.3554e-12,2.1032e-13. Six-mode cutoff guard; symmetry-zero entries use row scaling. Full metric completion is reported separately; this matching report alone has no source/flux conclusion.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_boundary_convergence.json`：Both actual original-author metric runs completed. Boundary controls6/5/rinf10000 ->7/6/rinf20000 at fixed xhor1e-4,WP50,goal24. Maximum author relative change at r2/r10/r30=6.83e-14/2.06e-14/3.34e-12, while local-vs-fine maximum remains8.191e-7. Hence the residual sampled local difference is not explained by these author boundary settings. This is not horizon-limit/high-L convergence or a final flux test.。

仍缺少：Isolated original physical kappa profiles/derivatives over the full actual source grid and boundary/high-ell convergence. Original operator/helper normalization and low-ell matching pass; complete spin0 contribution at three actual radii now agrees to1.16e-7, with a small r2 h_l-l- component increasing relative error. This is not yet the full environment source.

### Lorenz metric reconstruction and tetrad projection

实现：`src/lorenz_metric.py`, `src/lorenz_ghp.py`, `src/lorenz_tensor.py`, `src/paper_full_tetrad.py`, `src/environment_lorenz_mode.py`。

原文：DDKW eq:mp-s2,eq:mp-s1,eq:mp-s0,eq:h-projections; WDK2024 eq:hL and eq:H-DKW。

对照输入：all spin-sector radial jets and angular factors; a,r,theta,m,omega。输出：ten weighted tetrad components and BL metric,split sectors。

现有证据：
- `original_author_symbolic_code` — `outputs/lorenz_reference/Operators.nb`：Byte-identical WDK2024 source archive symbolic notebook; no local execution asserted。
- `independent_reimplementation` — `docs/environment_reproduction/full_tetrad_formula_audit.json`：Independent local derivation from published formulas。
- `internal_identity` — `docs/environment_reproduction/reconstruction_conditioning_audit.json`：High ell cancellation and perturbation sensitivity。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_schema.json`：Original 2026 a=.6,rp8,ell<=4 HDF5; several coefficients fail orbit-side continuity。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_metric_aligned_comparison.json`：Original-vs-local metric sample comparison; trace agrees, other differences need author-path/truncation diagnosis。
- `original_author_run` — `docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json`：This report validates the matching phase only; the completed actual-parameter full metric comparison is recorded separately below.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i6_h5_comparison.json`：Original Sam radiative reconstruction executed through all three spin sectors and the ten weighted metric projections at actual a=.8771530275949366,rp20,m1,L10,r=2,10,30. Protected ell1..3 component vector-relative maximum 8.191043e-07. Sector0/1/2 maxima 1.158013e-07,1.973502e-12,3.877619e-08. This is a completed full metric code path at three sampled radii and one m, not a full radial grid/high-ell/environmental source/flux verification. No scale/phase fitting; trace transformed by original Bmat0.。
- `original_author_input_propagation` — `docs/environment_reproduction/sam_actual_metric_i6_h5_source_propagation.json`：Saved actual original Sam/local metric coefficients, identically truncated to output spherical j1..3, reconstructed with spin-weighted spherical harmonics and converted to BL. m_g=-1 is conjugate of positive-m1 BL metric. Shared ThresholdCloud(alpha=.3),Hessian and scalar00 projection give J relative differences r2=2.4781e-9,r10=1.3254e-9,r30=1.8344e-7. Angular48->96 changes<=5.9e-15. Weighted-component cancellation max47.30 at r30 does not amplify this input error to percent scale. Neither result is fullJ or original author environmental-source code.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_i7_h6_comparison.json`：Completed original Sam full radiative metric code path, boundary-refined inford7/horord6/rinf20000, at actual a=.8771530275949366,rp20,m1,L10,r2/10/30 and protected spherical ell1..3. Each of ten component vectors compared over these three ell values: maximum relativeL2 8.191043e-07. Earlier i6/h5 author result changes by at most 3.336806e-12 under refinement. This is bounded full metric validation, not full-grid/high-ell/all-m/source/flux validation.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_metric_boundary_convergence.json`：Both actual original-author metric runs completed. Boundary controls6/5/rinf10000 ->7/6/rinf20000 at fixed xhor1e-4,WP50,goal24. Maximum author relative change at r2/r10/r30=6.83e-14/2.06e-14/3.34e-12, while local-vs-fine maximum remains8.191e-7. Hence the residual sampled local difference is not explained by these author boundary settings. This is not horizon-limit/high-L convergence or a final flux test.。
- `original_author_input_propagation` — `docs/environment_reproduction/sam_actual_metric_i7_h6_source_propagation.json`：Refined original author/local metric coefficients truncated identically to protected sphericalj1..3 and propagated through the SAME LOCAL ThresholdCloud/Hessian/scalar00 projection. Relative J differences r2/r10/r30=2.4781e-9/1.3254e-9/1.8344e-7;48/96 angular changes remain below6e-15. Weighted-component cancellation max47.30; no percent-scale amplification observed in these three truncated contractions. This is not original-author environmental source code and neither J is fullJ.。

仍缺少：Still missing original converged high-ell/full-source-grid metric data,additional m and original2025 environmental production settings. Both actual-parameter i6/h5 and refinedi7/h6 full metric code paths now pass ten component vectors at three radii for protectedell1..3 to8.2e-7. No remaining wait for the abandoned inford10 experiments.

### Source jumps and coupled matching

实现：`src/paper_all_component_matching.py`, `src/paper_sourced_matching.py`, `src/paper_jump_basis.py`, `src/paper_analytic_jumps.py`。

原文：DDKW eq:Jumps,eq:C0,eq:C1,eq:system。

对照输入：orbit data,selected spherical projections and spin2/trace jumps。输出：spin1/kappa jump constants; all metric continuity and derivative jumps。

现有证据：
- `independent_reimplementation` — `docs/environment_reproduction/all_component_matching_L8_q28_precise.json`：Independent matching solve in local code。
- `internal_identity` — `docs/environment_reproduction/full_gauge_closure_L8_precise.json`：Full closure check, no original jump arrays。
- `original_author_dataset` — `docs/environment_reproduction/author_hdf5_schema.json`：Direct raw orbit-side continuity diagnostic; no interpolation。
- `original_author_run` — `docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json`：Original Sam matching phase executed at a=.6,rp8,m1,L10,WP50. Four complex spin1/kappa jumps per ell1/2/3 agree with production at row-scaled1.98e-14,1.23e-12,1.15e-13. Modes ell4..10 are excluded as cutoff-edge diagnostics; not a full metric or flux comparison.。
- `original_author_run` — `docs/environment_reproduction/sam_actual_vs_production_jumps_ell1to3.json`：Original Sam matching executed at actual discrepancy a=.8771530275949366,rp20,m1,L10. Four spin1/kappa jump entries for ell1/2/3 agree with production at row-scaled errors 2.2993e-14,3.3554e-12,2.1032e-13. Six-mode cutoff guard; symmetry-zero entries use row scaling. Full metric completion is reported separately; this matching report alone has no source/flux conclusion.。
- `independent_reimplementation` — `docs/environment_reproduction/sam_vs_local_matching_L10_q32.json`：High-precision independent local L10/q32 matching assembled with original Sam lmlm supplementary equation compared to the executed original Sam a=.6,rp8,m1 result: ell1/2/3 relative2.33e-15,4.39e-16,4.01e-15. Same-matrix mixed-condition control4.76e-13,2.91e-14,4.03e-15. Cutoff-edge ell8/9/10 reach1.02e-8,1.70e-7,3.44e-6; only guarded ell1..3 validated. Absolute unused residual2.81 must not be reported as low-mode relative error.。

仍缺少：Full angular-cutoff convergence and high-ell matching remain. Original L10 matching passes guarded ell1..3 at both sample and actual parameters; first full metric code path has also completed at three actual radii, but this does not validate all high-ell/source-grid entries.

### Static sector and completion

实现：`src/environment_static_lorenz.py`, `src/lorenz_static_trace.py`, `src/lorenz_static_spin2.py`, `src/lorenz_static_gauge.py`。

原文：DDKW static eq:hm0,eq:Delta-c-def,eq:EG; Dyson completion discussion。

对照输入：a,rp,m=0,energy/angular momentum charge choices。输出：static metric and charges,regularity across boundaries。

现有证据：
- `internal_identity` — `docs/environment_reproduction/static_matching_convergence.json`：Local static refinement。
- `original_author_code` — `outputs/paper_metric_reference/srd24_KerrLorenzCirc/metric_reconstruction_calc_static.nb`：Author static notebook, not an execution record; special completion gauge convention described in README。

仍缺少：Raw author static completion data in same gauge at production parameters. Nonstatic gravitational m=1 produces dominant scalar m_s=0; these must not be confused.

### Covariant scalar source

实现：`src/environment_source.py`, `src/report_independent_scalar00_response.py`, `src/report_environment_forced_mode.py`, `src/source_provenance.py`。

原文：Dyson eq:operators2,eq:rhs_lead,eq:EOM_sourced: (Box-mu^2)phi11=h^{ab} nabla_a nabla_b phi10。

对照输入：normalized cloud,metric/Hessian,frequency mixing,angular quadrature,one-sided radial samples。输出：complex scalar source projected onto output spheroidal harmonics。

现有证据：
- `independent_reimplementation` — `docs/environment_reproduction/independent_scalar00_response_audit.json`：Alternative contraction/quadrature, but still shares our metric input。
- `internal_identity` — `docs/environment_reproduction/scalar_source_orbit_limits_comparison.json`：One-sided source limits。
- `original_author_code` — `docs/CONOR_ENVIRONMENT_SOURCE_SEARCH_20260915.md`：Complete fixed-commit readable-text search:2111files,2069texts,47.35MB,Git blob checked; no identifiable full massive cloud/Klein-Gordon environment source/flux pipeline found. Old notebooks contain real boson-tagged metric output paths, so this is not a claim of no environment-related artifacts. CompressedData/binary contents and unpublished local referenced files excluded.。
- `original_author_input_propagation` — `docs/environment_reproduction/sam_actual_metric_i6_h5_source_propagation.json`：Saved actual original Sam/local metric coefficients, identically truncated to output spherical j1..3, reconstructed with spin-weighted spherical harmonics and converted to BL. m_g=-1 is conjugate of positive-m1 BL metric. Shared ThresholdCloud(alpha=.3),Hessian and scalar00 projection give J relative differences r2=2.4781e-9,r10=1.3254e-9,r30=1.8344e-7. Angular48->96 changes<=5.9e-15. Weighted-component cancellation max47.30 at r30 does not amplify this input error to percent scale. Neither result is fullJ or original author environmental-source code.。
- `original_author_input_propagation` — `docs/environment_reproduction/sam_actual_metric_i7_h6_source_propagation.json`：Refined original author/local metric coefficients truncated identically to protected sphericalj1..3 and propagated through the SAME LOCAL ThresholdCloud/Hessian/scalar00 projection. Relative J differences r2/r10/r30=2.4781e-9/1.3254e-9/1.8344e-7;48/96 angular changes remain below6e-15. Weighted-component cancellation max47.30; no percent-scale amplification observed in these three truncated contractions. This is not original-author environmental source code and neither J is fullJ.。

仍缺少：Original author spline grids,derivative order,boundary extrapolation,source values and exact angular projection output; this is not supplied by public metric HDF5.

### Massive scalar radial Green solver

实现：`src/environment_radial.py`, `src/paper_leaver_radial.py`, `src/report_environment_forced_mode.py`。

原文：Dyson numerical appendix In/Up asymptotics,variation-of-parameters,Wronskian equations。

对照输入：mu,omega_s,l_s,m_s,angular eigenvalue,k,source grid,boundary/ODE controls。输出：massive In/Up radial functions,Green kernel,Z_H,Z_infinity。

现有证据：
- `independent_reimplementation` — `docs/environment_reproduction/leaver_radial_full_audit.json`：Our Leaver recurrence vs our ODE。
- `independent_reimplementation` — `docs/environment_reproduction/leaver_green_kernel_audit.json`：Independent kernel normalization。
- `internal_identity` — `docs/environment_reproduction/massive_radial_precision_audit.json`：21 fixed sources,endpoint/tolerance/Coulomb changes; not original code。
- `original_author_code` — `docs/CONOR_ENVIRONMENT_SOURCE_SEARCH_20260915.md`：Complete fixed-commit readable-text search:2111files,2069texts,47.35MB,Git blob checked; no identifiable full massive cloud/Klein-Gordon environment source/flux pipeline found. Old notebooks contain real boson-tagged metric output paths, so this is not a claim of no environment-related artifacts. CompressedData/binary contents and unpublished local referenced files excluded.。

仍缺少：Original author massive radial functions and normalization at threshold/superradiant regimes. Massless MST/GSN cannot be transferred unchanged to massive Klein-Gordon.

### Flux, mode sum and comparison target

实现：`src/environment_cloud.py`, `src/report_flux_coverage.py`, `src/report_full_reference_audit.py`, `src/extract_paper_flux_markers.py`。

原文：Dyson eq:enery_flux_ind,eq:Noether_ind,eq:scalar_flux,eq:energy_angmom。

对照输入：normalized Z_H,Z_infinity,omega_s,omega_c,m,cloud mass,mode coverage。输出：per-mode and summed energy/angular-momentum/Noether fluxes;paper plotted normalization。

现有证据：
- `plot_digitization` — `docs/environment_reproduction/full_reference_audit_20260915.json`：Quantitative comparison to PDF markers; not raw author arrays。
- `plot_digitization` — `docs/environment_reproduction/paper_v1_flux_markers.json`：Values extracted from original vector PDF figure。
- `original_author_document` — `outputs/environment_reference/Flux_Inf_Hor_TwoPanels_2.pdf`：Original author figure in arXiv source; graph coordinates are not original table。
- `original_author_code` — `docs/CONOR_ENVIRONMENT_SOURCE_SEARCH_20260915.md`：Complete fixed-commit readable-text search:2111files,2069texts,47.35MB,Git blob checked; no identifiable full massive cloud/Klein-Gordon environment source/flux pipeline found. Old notebooks contain real boson-tagged metric output paths, so this is not a claim of no environment-related artifacts. CompressedData/binary contents and unpublished local referenced files excluded.。

仍缺少：Original 2025 numeric per-mode fluxes,total flux arrays,mode cutoffs,error bars and exact settings for plotted points.

## 证据来源约束

- Original 2025 arXiv source contains the paper and figures, not cloud/source/radial/flux raw arrays or numerical solver.
- Original 2024 arXiv source contains amplitudes.dat and Operators.nb; these are real author reference artifacts.
- New srd24 and ConorDyson repository artifacts include real author source and HDF5, but current commits and 2026 fork sample are not demonstrated to be the 2025 production run.
- Original HDF5 trace has raw spheroidal basis; after explicit conversion to spherical, local trace agrees to <=9.75e-11 at six sampled points.
- HDF5 first-nine component discontinuities and inconsistent generator paths prevent using every number as ground truth.
- ConvolveSource.m from BHPTK/Teukolsky is a third-party upstream file, not Dyson adapted code; diagnose_*.py and src/paper_*.py are our implementations.
- Local metric/source/forced-mode caches and the patched MST bridge are not original author arrays/code.
- Original scalar/kappa operator helpers and original Sam matching have now actually run and passed bounded comparisons. The full original metric and environmental flux are separate milestones and must not be inferred from these module runs.
- Official BHPT Mathematica dependency code has been executed for11cases and matches AUTO radial R/Rprime to3.34e-14. These are official_dependency_run results, not the original environmental pipeline or a historical dependency lockfile.
- The first complete original Sam radiative metric execution at actual parameters has now finished. Ten protected ell1..3 components at three radii agree within8.2e-7; do not describe full metric accuracy as1e-14 or infer full environment source/flux validation.
- A full readable-text search of the Conor fixed commit found boson-tagged metric-output history but no identifiable complete massive environmental pipeline. This bounded absence of located code is not a claim that the authors do not possess it.
- The same bounded author/local metric inputs have been propagated through a SHARED LOCAL cloud/Hessian source. J differences<=1.84e-7 at three radii rule out percent-scale amplification for that truncated diagnostic only; they do not validate original author cloud/source code or the complete source.
- Both i6/h5 and i7/h6 original full metric sample runs completed; boundary refinement changes author values by<=3.34e-12 but local-vs-author error remains<=8.2e-7. Inford10 symbolic-expansion experiments were abandoned as complete-run candidates, with their completed matching outputs retained; they are not waiting tasks.

完整分类、来源提交、arXiv包成员与哈希见同名JSON。共享源传播脚本为`src/report_author_metric_source_propagation.py`；粗/细结果分别为`sam_actual_metric_i6_h5_source_propagation.json`和`sam_actual_metric_i7_h6_source_propagation.json`。
