# Kerr 点粒子数值 solver 与数据迁移计划

- 来源任务：T007 / `TASK_num.md` 任务8
- 状态：`completed`；owner 批准的 M0--M7 迁移已执行
- Owner 决定：W02 已结束，当前理论、符号验证和数值结果可信，可以晋升
- 既有迁移：[T006 理论与符号证据迁移记录](../findings/kerr_point_particle_evolution/migration_plan.md)
- 历史来源：原 W02 数值 workspace 已迁出仓库；历史计划仍在 [W02 plan tree](../workspace/plan/W02_point_particle_evolution/)
- 历史计划：[W02 plan tree](../workspace/plan/W02_point_particle_evolution/)

## 1. 目标与迁移边界

本次只处理 T006 之后完成的 W02.3 数值部分。连续理论已经位于 `findings/kerr_point_particle_evolution/`，Mathematica checks/results 已经位于 `code/kerr_point_particle_verification/` 与 `data/kerr_point_particle_verification/`；这些正式证据不重复移动，也不作为“历史验证数据”删除。

迁移遵循：

1. 已审阅代码和必要测试优先直接移动，不在搬运过程中重写算法；
2. 正式代码、科学结论、最终数据和历史计划分别进入 `code/`、`findings/`、`data/`、`workspace/plan/`；
3. 结论提取并完成 hash/loader 检查后，删除 smoke、checkpoint、cache、重复 run 和过时审阅材料；
4. SXS:BBH:0305 的 $(2,2)$、$(4,4)$ medium 结果作为两个不可删除的主数据 bundle；
5. 迁移后再解决两个已知 workflow 缺口，不能用迁移掩盖尚未实现的接口。
6. 迁移后的新计算直接写入 `data/kerr_point_particle_evolution/` 的受控 staging 路径，不再以 workspace run 作为正常中转。

Owner 对 W02 的接受允许正式迁移，但不应篡改已有 artifact 的事实性 metadata。例如现有 C2 文件中的 A1 sign/normalization `conditional-open` 应忠实保留；位置变化本身不把未独立检查的 convention 改写成已验证。

## 2. 当前内容分类

### 已在 T006 晋升

以下内容保持原位：

- `findings/kerr_point_particle_evolution/theory/`：主推导、约定与全部理论附录；
- `code/kerr_point_particle_verification/checks/`：25 个 Wolfram Language checks；
- `data/kerr_point_particle_verification/results/`：机器 residual；
- `findings/kerr_point_particle_evolution/evidence/theory_status_audit.md`。

这些是正式理论的可复现证据，不是可删除的临时 validation output。

### 本次直接移动

| Workspace 来源 | 正式目标 | 处理 |
|---|---|---|
| `numerical/src/sminus2_point_particle/` | `code/kerr_point_particle_solver/src/sminus2_point_particle/` | 整体移动；保持 package 名与公式实现 |
| `numerical/pyproject.toml` | `code/kerr_point_particle_solver/pyproject.toml` | 直接移动 |
| 核心 `scripts/` | `code/kerr_point_particle_solver/scripts/` | 直接移动；清除绝对 workspace 路径 |
| `tests/test_formula_contracts.py`、`test_c1_c2.py`、`test_strain.py` | `code/kerr_point_particle_solver/tests/` | 与实现一起移动 |
| `tests/fixtures/chi0p8_l6_flux.json` | `code/kerr_point_particle_solver/tests/fixtures/` | 作为回归 fixture 保留 |
| `validation/formula_code_map.md`、`numerical_contract.md` | `code/kerr_point_particle_solver/docs/` | 直接移动后更新正式路径 |
| `validation/c2_implementation.md` | `code/kerr_point_particle_solver/docs/strain.md` | 保留公式、实现合同、主审阅与边界 |
| 当前 `README.md` | `code/kerr_point_particle_solver/README.md` | 以现有内容为基础更新安装、运行和正式状态 |

`benchmark_rhs.py` 随 solver 保留；基础 $\psi_4/H$ 绘图进入公共工具。历史 `compare_mode_waveforms.py` 不满足 W03 的 phase-only 科学合同，结论迁入 findings 后删除。`compute_isco_flux.wl`、`write_flux_request.py` 和 `build_initial_data.py` 作为新自旋首次生成初值的正式入口保留。

### 提取后删除

以下内容不晋升为正式 artifact：

- `.DS_Store`、`__pycache__/`、`.pytest_cache/`；
- 所有 smoke run 及其 checkpoints；
- $\chi=0.8$ 跨代码 comparison run 中的 checkpoints、重复 waveform 和绘图中间件；
- SXS:BBH:0305 两个已完成 run 的 rolling checkpoints；
- `psi4_pipeline_status_review.md`、旧 `c1_c2_implementation_status.md`、`review_status.md` 等已被最终结论取代的过程性审阅；
- 重复的 smoke/compare 配置；
- 临时 profiling 输出和可以由正式数据重画的重复图片。

跨代码对比的结论、关键数值与接受理由写入 findings 后，原 comparison 数据整体删除；不把历史验证 run 当作第三套正式波形保存。性能优化只在 solver README/docs 保留最终 benchmark 摘要，过程性记录删除。

### 必须留在 workspace

`workspace/plan/W02_point_particle_evolution/` 全树保留，作为 W02 的历史计划和 owner 决策记录。迁移结束后，`workspace/W02_point_particle_evolution/README.md` 只保留简洁的完成状态与正式目录链接；不保留第二份 solver 源码或正式数据。

## 3. 建议正式文件树

~~~text
findings/
└── kerr_point_particle_evolution/
    ├── README.md
    ├── numerical_evolution.md
    └── evidence/
        └── cross_code_low_resolution.md

code/
├── kerr_point_particle_verification/       # T006 已有，保持
├── kerr_point_particle_solver_migration.md # 本计划；执行后更新为迁移记录
├── kerr_waveform_tools/
│   ├── AGENTS.md
│   ├── README.md
│   ├── pyproject.toml
│   ├── src/kerr_waveform_tools/
│   │   ├── contracts.py
│   │   ├── background.py
│   │   ├── transition.py
│   │   ├── worldline.py
│   │   ├── harmonics.py
│   │   ├── artifacts.py
│   │   ├── waveform_io.py
│   │   └── strain.py
│   ├── scripts/
│   │   ├── compute_isco_flux.wl
│   │   ├── write_flux_request.py
│   │   ├── build_initial_data.py
│   │   ├── plot_mode_real_abs.py
│   │   └── plot_strain.py
│   └── tests/
└── kerr_point_particle_solver/
    ├── AGENTS.md
    ├── README.md
    ├── pyproject.toml
    ├── src/sminus2_point_particle/
    ├── scripts/
    ├── tests/
    │   └── fixtures/
    ├── configs/
    │   ├── smoke/
    │   └── examples/
    └── docs/
        ├── formula_code_map.md
        ├── numerical_contract.md
        ├── performance.md
        └── strain.md

data/
├── kerr_point_particle_verification/       # T006 已有，保持
└── kerr_point_particle_evolution/
    ├── runs/
    │   ├── .staging/<run_id>/              # incomplete / restartable
    │   └── <run_id>/                       # complete，仍可 awaiting_review
    └── sxs0305_medium/
        ├── manifest.json
        ├── input/
        │   ├── sxs_metadata.json
        │   ├── sxs_h22_h44.npz
        │   ├── bhpt_l6_isco_flux.json
        │   └── plunge_initial_data.json
        ├── mode_22/
        │   ├── run_metadata.json
        │   ├── psi4_l2_m2.npz
        │   ├── H_l2_m2.npz
        │   ├── psi4_real_abs.png
        │   └── strain_real_abs.png
        └── mode_44/
            ├── run_metadata.json
            ├── psi4_l4_m4.npz
            ├── H_l4_m4.npz
            ├── psi4_real_abs.png
            └── strain_real_abs.png

workspace/
├── plan/W02_point_particle_evolution/      # 历史计划完整保留
└── W02_point_particle_evolution/README.md  # 完成导航
~~~

### 新计算的输出路径合同

Repository wrapper 的默认 output root 固定为

~~~text
data/kerr_point_particle_evolution/runs/
~~~

新 case 先创建

~~~text
data/kerr_point_particle_evolution/runs/.staging/<run_id>/
~~~

并在其中保存 `status=incomplete` 的 metadata、rolling checkpoints 和分块数据。只有 PDE 演化、请求的 C2 转换、两张基础图、hash 和 loader 检查全部完成后，才把整个目录在 `data/` 内原子定稿为 `runs/<run_id>/`。定稿状态 `complete` 只表示计算完整；是否被 owner 接受仍由独立 `review_status` 表示，不能由目录位置推断。

可复用 Python API 仍允许调用者显式提供其他 output root，但仓库内正式 CLI 和 production config 默认直接使用上述 `data/` 路径。代码、README、示例和 metadata 不再把 `workspace/W02_point_particle_evolution/numerical/runs/` 写成默认输出位置。

每个默认正式 fixed-$(\ell,m)$ case 应同时请求并保留

$$
\psi_{4,\ell m}(T),\qquad H_{\ell m}(T)=r h_{\ell m}(T),
$$

并在同一 bundle 自动生成：

- `psi4_real_abs.png`：$\operatorname{Re}\psi_{4,\ell m}$ 与 $|\psi_{4,\ell m}|$；
- `strain_real_abs.png`：$\operatorname{Re}H_{\ell m}$ 与 $|H_{\ell m}|$。

两图都覆盖完整实际时间序列，标出 $T_{\rm stable}$、light-ring reference、source-off 和 FFI taper/推荐区间。图片是标准交付，不替代对应 complex NPZ。若显式关闭某一数据产品，只允许用于 test/diagnostic run；production default 同时保存 $\psi_4$、strain 及两张图。

## 4. SXS:BBH:0305 正式数据

### 背景与输入

SXS:BBH:0305 使用本地 `SXS:BBH:0305v2.0/Lev6`、N2 数据：

$$
\frac{M_f}{M}=0.952032939704,
\qquad
\chi_f=0.6920851868180025.
$$

精确自旋的 $\ell_{\max}=6$ ISCO flux 已由 Mathematica MCP 的 `KerrGeodesics`/`Teukolsky` machine-real 路径计算：

$$
\frac{\dot E_{\rm ISCO}}{\mu^2}
=0.009328435794991877,
\qquad
\frac{\dot E_{\ell=6}}{\dot E_{\ell\le6}}
=0.011484701037361529.
$$

由此得到

$$
r_0=3.4157225309091355M,\qquad
E_p=0.8975004479262404,\qquad
L_{z,p}=2.6006748443128864M.
$$

Formal manifest 必须保存 SXS version/Lev/N2、SXS cache hash、$M_f,\chi_f$、flux 值、Wolfram 14.2.1/MacOSX-ARM64、初值 key 和生成脚本路径。0305 不属于现有 W01 22-case manifest，因此建立独立 comparison bundle，不静默修改 W01 catalog。

### 两个不可删除的主 bundle

| Bundle | 当前原始数据 | 固定设置 |
|---|---|---|
| $(2,2)$ | `runs/sxs0305_medium_m2_20260903/psi4_l2_m2.npz` 与 `metadata.json` | $m=\ell=2$ |
| $(4,4)$ | `runs/sxs0305_medium_m4_20260903/psi4_l4_m4.npz` 与 `metadata.json` | $m=\ell=4$ |

共同设置为

$$
(N_R,N_y)=(512,129),\quad
\Delta T_{\rm out}=0.1M,\quad
\sigma_R=4\Delta R,\quad
\sigma_y=4\Delta y,\quad
\tau_{\rm on}=20M,
$$

$$
T_1=358.7M,\qquad
N_T=3587,\qquad
\Delta T_{\rm post}=120M.
$$

每个 bundle 保留原始 complex $\psi_{4,\ell m}$ 作为权威 solver 输出，再由正式 C2 生成完整 $H_{\ell m}=rh_{\ell m}$；strain 不能替换或覆盖 $\psi_4$。迁移前后分别计算 SHA256，并用 loader 检查 mode、background、sample count、uniform $T$、finite values、status、code/config/trajectory hashes。

现有 run 的 `code_hash=27b96348...` 是生成时的整包 hash。M0 实测当前 package hash 与它完全一致，因此迁移 manifest 需记录：

- 原始 generating package hash；
- 晋升代码的 release hash；
- 两者在 M0 的一致性结论；
- 迁移后公共代码抽取与 workflow 修复产生的新 release hash，不回写原 run metadata。

## 5. 可跨 W03--W06 复用的公共代码

### 分层原则

后续任务不应为了调用 ISCO、worldline、waveform loader 或 FFI 而依赖名为 `sminus2_point_particle` 的完整时域 solver。公共代码应满足：

$$
\text{kerr\_waveform\_tools}
\;\not\longrightarrow\;
\text{kerr\_point\_particle\_solver},
$$

而 solver 可以单向依赖公共 package。公共层只接收明确的物理对象、数组和 provenance，不接收 `SpatialGrid`、`SpatialRHS` 或 run-global 单例。

为了遵守“升级内容优先直接移动”，第一步仍把现有 package 原样移动到正式 solver；随后在单独阶段做等价抽取。抽取完成后只有公共 package 保存真实实现，`sminus2_point_particle` 在一个兼容周期内只 re-export 旧 API，不保留第二份函数副本。

### 现有代码的公共 API 候选

| 当前对象 | 公共目标 | 主要使用者 | 拆分要求 |
|---|---|---|---|
| `errors.ContractError` | `contracts.py` | W03--W06、solver | 统一输入/schema 错误；solver-only `OpenBoundaryError` 留在 solver |
| `geometry.delta/sigma/horizons`、compactification、`isco_quantities`、`RadialState/radial_state` | `background.py` | W05 近视界几何、W06 phase/radial inputs、所有新轨迹 | 保留显式 $(M,a,L)$；不绑定 W02 grid |
| `trajectory.prograde_light_ring_radius`、`tortoise_radius` | `background.py` | W03 的 $u_{\rm LR}$、W05/W06 | 与 `T|_{\mathscr I^+}=u` 的约定链接到 findings |
| `TransitionManifest`、`PlungeInitialData`、`FluxRecord`、Ori--Thorne generator、versioned cache | `transition.py` | 后续每个 Kerr case、W05/W06 worldline | 去掉对 `EvolutionConfig` 的反向依赖，改收显式 background/transition spec |
| `Jet3`、`StageSample`、`trajectory_rates`、`DOP853Trajectory`、`integrate_trajectory` | `worldline.py` | W03 LR 时间、W05 直接积分、W06 saddle phase | 保留唯一 $(R_p,\Phi_p)$ 与事件定义；不导入 solver source/grid |
| `spin_weighted_spherical_harmonic`、`projection_row`、`project_mode` | `harmonics.py` | W03/W04 mode convention、W05 spherical projection checks | `scri_slice` 与 grid-specific Simpson assembly 留在 solver；不声称提供 spheroidal source |
| `ArtifactStatus`、count/hash、atomic JSON/NPZ primitives | `artifacts.py` | W03--W06 的 cache、tables、waveforms | 从当前 `io.py` 中拆出，不依赖 `EvolutionState` |
| complex mode/H schemas 与只读 loaders | `waveform_io.py` | W03 unfiltered、W04 filter、W05/W06 comparisons | Solver checkpoint/full-field I/O 留在 solver；v1/v2 compatibility 在公共 loader |
| FFI、taper、cutoff/sampling diagnostics | `strain.py` | W03 的共同 strain、W04 filter 输入 | 保持 $H_{\ell m}=rh_{\ell m}$、完整 transform 与 conditional convention metadata |
| `compute_isco_flux.wl`、`write_flux_request.py`、`build_initial_data.py` | 公共 `scripts/` | 所有新 Kerr backgrounds | 作为 transition/flux 数据生成入口，不藏在 solver 示例中 |
| 基础 complex mode 绘图逻辑 | `plotting.py`、`plot_mode_real_abs.py`、`plot_strain.py` | W02 数据 bundle、W03--W06 快速审阅 | 参数化 field label/events；不硬编码 $\chi$、$(\ell,m)$ 或窗口 |

`geometry.py`、`initial_data.py` 和 `trajectory.py` 是同一依赖链，不能只移动 `isco_quantities` 而留下两套 background/transition 定义。`io.py` 必须真正拆分：generic artifact/waveform I/O 进入公共层，`Checkpoint`、`RollingCheckpoint`、`FullFieldMemmap` 等 solver state 留在 W02 solver。

### 后续任务的复用关系

| 工作 | 直接复用 | 不在本次提前实现 |
|---|---|---|
| W03 | waveform/H loaders、background mass metadata、$r_{\rm LR}$、$r_*$、trajectory events、harmonic convention | 实际质量/时间对齐、phase-only fit、mismatch scan、QNM table |
| W04 | waveform artifacts、uniform-series/taper 基础、atomic outputs | rational QNM filter、filter lists 与 control injections |
| W05 | background、ISCO/transition、完整 worldline、四速度/time jets、spherical-harmonic projection checks、cache primitives | spheroidal angular source、`pybhpt` radial cache、近视界展开与直接频率积分 |
| W06 | background/worldline、artifact cache、参数与环境 provenance | complex-saddle solver、WolframClient radial requests、Hessian/contour logic |

QNM、rational filter、`pybhpt` radial cache 和 SPA 尚不存在于当前已审阅 W02 code，不能为了填满公共库而预先创建空 abstraction。它们在各自任务形成第一份可执行实现后，再依据真实的第二个消费者决定是否进入公共 package。

### 不应抽到公共层

- `field.py`、`source.py`：完整 W02 fixed-$m$ hyperboloidal Teukolsky operator/source；
- `grid.py`、`operators.py`、`evolution.py`：W02 method-of-lines、boundary rows 与 RK4；
- `workflow.py`、`EvolutionConfig`、`OutputRequest`：W02 solver assembly、运行和 direct-to-data routing；
- `benchmark_rhs.py`：W02 性能工具；
- 当前 `compare_mode_waveforms.py`：允许一般 complex scale，不满足 W03 的 phase-only 科学合同，只作为旧 solver validation helper 保留或删除，不能冒充 W03 公共 comparison API。

### 抽取验收

- `kerr_waveform_tools` 可以独立安装和运行 tests，不导入 `sminus2_point_particle`；
- solver 依赖公共 package 的正式版本，不使用相对路径注入；
- 旧 `sminus2_point_particle` 公共名字只通过显式 compatibility re-export 存活，并有 deprecation 说明；
- 抽取前后 ISCO、transition、trajectory events、harmonics、FFI 和 loader fixtures 保持相同结果；
- 公式权威仍是 findings；公共化不扩大已验证的物理适用范围。

## 6. 迁移后必须解决的两个实现缺口

直接移动后识别出的两个 workflow 缺口均已关闭；本节保留问题与修复合同。

### 6.1 `save_strain_lm` 必须兑现

当前 `workflow.py` 只用 `save_strain_lm=True` 决定是否形成 $\psi_{4,\ell m}$，不会调用 FFI 或写出 $H_{\ell m}$。正式接口应满足：

| `save_psi4_lm` | `save_strain_lm` | 最终产品 |
|---|---|---|
| false | false | 不生成 mode/strain |
| true | false | $\psi_{4,\ell m}$ |
| false | true | $H_{\ell m}$；$\psi_4$ 仅作可恢复的内部 staging |
| true | true | 同时保存 $\psi_{4,\ell m}$ 与 $H_{\ell m}$ |

FFI 只在完整 accepted run flush 后执行，并使用完整 mode series：

$$
T_{\rm FFT}=T_{\min},\qquad
T_a=T_{\rm stable},\qquad
T_b=T_1.
$$

若只请求 strain，为支持 restart，可把 $\psi_4$ 写入 run 内明确标记的 staging artifact；只有 $H_{\ell m}$ 原子写入并由 loader 复核后才删除 staging。中断或转换失败时保留可恢复数据，不写 `complete` strain marker。

需要新增四种 output-request 组合、test/full、restart、FFI failure 和多 $\ell$ 的测试。该行为改变 run output contract，新生成数据应升级到 `c1-run-v2/c1-mode-v2`；旧 SXS0305 `v1` 数据只读兼容，不重写。

### 6.2 `run_case.py` 必须消费 versioned initial-data cache

当前入口每次从 flux JSON 调用 `generate_ori_thorne_initial_data()`，没有查询 `JsonInitialDataCache`。正式流程应为

$$
\mathcal K_{\rm ID}
\longrightarrow
\begin{cases}
\text{validated cache hit}, & \text{直接读取},\\
\text{cache miss}, & \text{由 versioned flux 生成、验证并原子写入}.
\end{cases}
$$

CLI 提供明确的 cache directory；run metadata 记录 key、hit/miss、record path/hash、schema/generator/flux versions。Cache hit 仍需重算关键一致性，不因文件存在而信任。测试覆盖首次 miss、再次 hit、key mismatch、tampered record、不同 $m$ 共用同一背景初值，以及只读 cache 失败。

SXS0305 的 $m=2$ 与 $m=4$ 必须解析为同一 initial-data key 和相同 $(r_0,E_p,L_{z,p})$，但 run config/hash 保持各自独立。

## 7. 分阶段执行

### M0. Freeze 与 inventory

1. 确认两个 SXS0305 run 为 `complete`，并由 loader 读取；
2. 生成迁移 inventory：相对路径、大小、SHA256、artifact/schema/status；
3. 运行当前完整 tests，记录基线；
4. 记录旧 run hash 与当前 release hash 的差异边界；
5. Owner 审阅本计划后才进入 M1。

### M1. 原样晋升

建立 `code/kerr_point_particle_solver/`，直接移动 package、scripts、tests、fixtures 和核心 docs。只允许机械更新 imports、相对路径、README 状态与安装入口；不在这一阶段修改数值方法。

### M2. 公共代码等价抽取

按第 5 节依赖边界建立 `code/kerr_waveform_tools/`。先拆 contracts/background/transition/worldline，再拆 artifacts/waveform I/O/strain/harmonics；每一步让 solver 改为正式 import，并由原 fixture 做前后比较。不能复制后再让两份实现同时存活。

### M3. Workflow 缺口修复

先完成 initial-data cache hit/miss，再完成 `save_strain_lm` 四态路由和 direct-to-data staging/finalize router。每项单独提交实现与测试；不改变 PDE、source、grid、RK4、Gaussian 或 C2 数学公式。启用 v2 schema，并保留 v1 read-only loader。Production default 同时启用 $\psi_4$、strain 和两张基础图；plot failure 不得留下伪 `complete` bundle。

### M4. SXS0305 数据 bundle

1. 把两个 medium $\psi_4$ 和原 run metadata 直接移动到 data；
2. 保存 SXS h22/h44 输入及其小型 metadata；
3. 保存精确 flux 和 plunge initial-data record；
4. 用修复后的正式 C2 分别生成 $H_{22}$、$H_{44}$；
5. 为两个 mode 分别生成完整 `psi4_real_abs.png` 与 `strain_real_abs.png`；
6. 写 `manifest.json`，逐项记录 source/code/config/trajectory/input/output hashes；
7. 比较所有直接移动数组的迁移前后 hash。

### M5. Findings 与导航

新增 `numerical_evolution.md`，只写已接受的物理问题、方程对象、默认离散、跨代码结论、0305 数据和 C2 边界；不复制 implementation debug history。更新既有 findings README、PROJECT task/file tree 和相关链接。

### M6. 清理

在 M0--M5 全部通过后，删除第 2 节列出的缓存、smoke、checkpoint、重复 comparison 数据、过时 validation notes 和冗余 configs。清理后再次搜索旧路径、竞争性 README、绝对 workspace 引用和 stale status。

### M7. 任务收尾

运行正式 package tests、两个短 CLI、cache miss/hit、四种输出路由和 data loaders；记录最终 inventory。随后更新 `work_log.md`。只有 owner 明确要求时才执行 Git commit，并关联 T007。

## 8. AGENTS 与 FutureSKILL proposal

以下五项是仓库组织的建议方向，不是固定模板。实际落地时应服从具体任务、已有结构和更清楚的表达；当例外能提高科学可审阅性、可复现性或工程可靠性时，可以采用例外并简要说明原因。

| Proposal | 候选位置 | 建议倾向 | FutureSKILL 中的通用表达 |
|---|---|---|---|
| 1. 晋升与修复适当分开 | `code/AGENTS.md` | 对已经审阅的实现，通常优先保持内容等价地移动；若迁移同时需要行为修复，倾向于把修复及其测试单独说明，使算法变化容易辨认 | Promotion 通常保持 reviewed artifact 等价；必要的行为变化应保持可见并重新审阅 |
| 2. Direct-to-data lifecycle | `data/kerr_point_particle_evolution/AGENTS.md`，必要时在 `data/AGENTS.md` 提供入口 | 对高成本正式计算，可直接使用 data 下的 staging/finalize 流程；通常区分计算完成与 owner review，并随最终 bundle 保留足以复现和审阅的 complex 数据、metadata、provenance 与基础图。Checkpoint 是否删除由恢复需求决定 | 高成本任务可采用 direct-to-data staging；目录位置通常不自动代表 approved |
| 3. 历史计划与正式产物分工 | `workspace/AGENTS.md` | 任务完成后，倾向于保留仍有解释价值的 plan 和简洁导航，并减少与正式代码、数据竞争的重复副本；若 workspace 副本仍有独立用途，可以明确标注后保留 | Plan 可以作为决策记录保留；正式 artifact 倾向于有清楚的权威入口 |
| 4. Solver 本地边界 | `code/kerr_point_particle_solver/AGENTS.md` | Solver 通常从 findings 读取公式基准，把运行数据放在 data；可能改变科学结果的数值方法倾向于先回到 workspace 审阅。I/O/schema 演进时，按历史数据价值决定是否提供兼容读取 | 项目特定，通常不提升到通用 FutureSKILL |
| 5. 公共 Kerr 工具边界 | `code/kerr_waveform_tools/AGENTS.md` | 公共层倾向于避免依赖完整 W02 solver，并让单位、convention 和 provenance 在接口中可见；公共 abstraction 通常由真实的跨任务使用需求驱动，不必提前为尚未实现的算法预设结构 | 公共代码应更多由真实消费者推动，避免从单一 workflow 过度抽象 |

实际执行时先检查现有规则并只更新受影响的位置。FutureSKILL 可吸收 proposal 1--3 与 5 的通用倾向，但不必复制 SXS、Teukolsky、固定目录名或本项目的具体例外。

## 9. 可复用 prompts

以下 prompts 从 T007 任务1--8的完整过程提取，分别覆盖数值规划、实现执行和任务结束后的迁移。执行迁移时，`PROMPTS.md` 最多保留这三个通用版本；不保存只对某个小节有效的长对话原文。

### Prompt 1：从物理流程建立数值计划

> 依据已经固定的连续方程和物理 endpoints，建立一份物理学家可读、又能指导代码建设的数值计划。先用简洁流程图说明从背景/轨道参数，经网格、初值、trajectory、source、method-of-lines 和时间推进，到 SCRI+ mode 与后处理输出的主链；再按同一顺序写实施步骤。主计划强调物理对象、公式、模块关系和数据流，有限差分边界、CFL、RK 内部步长与输出 cadence、startup、I/O/restart 等细节按需要放入链接附录。保持已有可靠框架，不为填模板增加内容；未确定项和 owner review 点明确标记。

### Prompt 2：按计划实施、审计并运行数值代码

> 以已审阅 numerical plan 为合同，按物理流程逐步实现代码。先固定公式到代码的映射和默认参数，再实现 trajectory、grid/boundary、source/RHS、RK 时间推进、输出/restart 和后处理；区分 $\Delta T_{\rm RK}$ 与 $\Delta T_{\rm out}$，每个 RK stage 在正确时间采样动态 source。Coding agent 完成后，由指定审计角色检查 plan-to-code 一致性、数值危险、I/O/status/provenance 和失败路径；只在出现会改变物理结论或计划结构的问题时暂停请求 owner。长计算先做 smoke，提供进度，使用 profiling 优化等价代码而不放宽数值方法，并用完成标记、loader 和必要的独立结果比较判断交付。

### Prompt 3：结束任务并迁移可信科研产物

> Owner 接受当前任务后，先盘点 workspace、现有正式目录和后续任务的真实消费者，给出直接移动、公共代码抽取、findings 提炼、最终数据 bundle、历史计划保留和删除候选的完整映射。已审阅内容倾向于先等价移动；行为修复与公共 API 抽取另立阶段并重新测试。高成本新计算可采用 direct-to-data staging/finalize，保留原始与派生 complex 数据、基础图、配置、事件、环境和 hashes，同时区分 complete 与 owner review。迁移前后核对 hash/loader/schema/链接，保留有价值的计划，确认正式产物完整后再删除 checkpoint、smoke、cache、重复验证数据和过时摘要；最后同步 AGENTS、PROJECT、work log、FutureSKILL 和不超过任务要求数量的通用 prompts。

## 10. 验收条件

迁移完成需同时满足：

- 正式 solver 可从独立 checkout 安装、运行和测试，不依赖 workspace 绝对路径；
- `kerr_waveform_tools` 可独立安装，W03--W06 所需的现有公共对象不再依赖完整 W02 solver；
- 仓库正式 CLI 默认直接写入 `data/kerr_point_particle_evolution/runs/.staging/`，成功后在 data 内原子定稿；
- `save_strain_lm` 四态输出与 initial-data cache hit/miss 已实际实现；
- v1 SXS0305 数据保持原 metadata/hash，可由正式 loader 读取；新 run 使用 v2 schema；
- SXS0305 的 $(2,2)$、$(4,4)$ 两个 bundle 同时含原始 $\psi_4$、完整 C2 strain、两张标准波形图、run/input provenance 和 manifest；
- 理论、MMA、solver、findings 与 data 之间只有一套权威副本；
- 历史计划完整留在 workspace，过程性验证数据和缓存已删除；
- 不超过五项 AGENTS proposal 与三个 prompts 已按 owner 决定落地；
- PROJECT、README、work log、FutureSKILL、链接和文件树与最终结构一致；
- Git diff 不含意外公式、参数、边界条件或数值方法变化。

## 11. 执行结果

M0--M7 已按 owner 决定完成：

- `code/kerr_point_particle_solver/` 保存 35 个 solver、test、config 与文档文件；`code/kerr_waveform_tools/` 保存 20 个公共文件，公共包不依赖完整 solver；
- `save_strain_lm` 四态输出、`c1/c2` v1/v2 兼容读取、initial-data cache miss/hit 和 direct-to-data staging/finalize 均已实现；
- 正式 SXS0305 bundle 保存 20 个文件，其中 manifest 跟踪 16 个输入/结果文件；两个原始 $\psi_4$ 与 run metadata 的迁移前后 SHA256 不变；
- solver tests 为 `58 passed`，公共包 tests 为 `4 passed`；两个真实短 CLI 分别得到 `miss-generated` 与 `hit`；manifest、四个 mode/strain loaders 和源码 hash 均通过；
- SXS0305 $(2,2)/(4,4)$ medium 的原始 $\psi_4$、完整 FFI $H$ 与四张标准图进入 `data/`，跨代码结论进入 findings；
- 历史 plan 保留；重复 smoke/comparison/checkpoint 与过程性审阅已从仓库移出。安全策略不允许直接不可恢复删除整个旧数值目录，因此迁移时将其暂存于 `/private/tmp/t007_workspace_numerical_pre_cleanup_20260904`；
- AGENTS、FutureSKILL、PROJECT、README 和三个 T007 reusable prompts 已同步。未执行 Git commit；是否提交仍由 owner 决定。
