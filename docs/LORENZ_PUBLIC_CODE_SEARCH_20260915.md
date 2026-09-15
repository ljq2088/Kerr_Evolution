# Kerr–Lorenz 原作者公开代码检索与模块核验（2026-09-15）

本轮找到了此前未使用的公开原作者代码。原先“只能向作者索取、没有公开 Kerr 度规重构实现”的判断已经过时。现已固定下载 Sam Dolan 的实现和 Conor Dyson 的后续分支；但是，**后续分支含有不同版本的低阶截断和调试 loader，不能把其任意文件视为 2025 环境论文的原始真值。**

这份报告只负责公开源码、变量约定和独立执行的基础模块核验。公开 HDF5 与本地全十个度规分量的实际对比由同轮单独诊断负责；本报告不重复运行或预先宣称其结果。

## 来源和可用范围

| 来源 | 固定提交 | 获取内容与用途 |
|---|---|---|
| [Sam Dolan / KerrLorenzCirc](https://github.com/srd24/KerrLorenzCirc/tree/5c1b42793ff893fd0c65a5b48ceecc02d34c5e35) | `5c1b42793ff893fd0c65a5b48ceecc02d34c5e35` | 62 个文件；原始 Mathematica notebook、生成的 radiative `.m`、静态模式和 completion |
| [Conor Dyson / KerrLorenzMSF](https://github.com/ConorDyson/KerrLorenzMSF/tree/6b469f0835e42d3fa531b4997d7a912257119f60) | `6b469f0835e42d3fa531b4997d7a912257119f60` | 2111 个文件；后续模块化 `.wl`、MST/NGrid 接口、HDF5 数组、实验/旧目录 |

GitHub 元数据记载 Sam 仓库创建于 2025-04-07，Conor 分支创建于 2026-01-15；下载提交分别最后推送于 2025-04-07 和 2026-06-07。这些日期是仓库元数据，不能证明算法第一次编写或 HDF5 生成日期。

Sam 的 README 明确说明其实现 [2023/2024 圆赤道轨道 Lorenz 重构工作](https://arxiv.org/abs/2306.16459) 与 [2021/2022 Kerr Lorenz 重构工作](https://arxiv.org/abs/2108.06344)，并致谢 Conor Dyson 参与测试。Conor 仓库没有根目录 README，因而不能将它的现有 HDF5 直接归属于 2025 环境论文。Toolkit 的 [h1Lorenz](https://github.com/BlackHolePerturbationToolkit/h1Lorenz) 项目限于 Schwarzschild，不能代替这里的 Kerr 参考。

完整下载和逐文件 SHA-256 清单在 `outputs/paper_metric_reference/manifest.json`。精简且适合纳入版本控制的来源、关键文件哈希、源码行和实际执行结果在 `docs/lorenz_public_code_manifest_20260915.json`。原始下载未修改。

## 实际可比模块

Sam `metric_reconstruction_calc_radiative.m` 包含角向球谐展开、解析曲率与迹跳跃、spin-1/κ 的源匹配线性系统、h/κ 的径向 ODE 与边界级数、十个 tetrad 分量重构及球谐投影。静态另有三个构建 notebook 与 completion 文件。它已覆盖此前最缺少的 **低阶规范势、源匹配、κ particular 与齐次补偿** 技术。

| 作者模块 | 本地对应 | 必须对齐的内容 |
|---|---|---|
| 自旋椭球谐、球谐投影与 Γ 矩阵 | `paper_precise_angular.py`、`environment_angular_diagnostic.py` | SWSH 相位、归一化、Chandrasekhar 的负自旋本征值约定 |
| 解析曲率跳跃、迹跳跃 | `paper_analytic_jumps.py`、`paper_sourced_matching.py` | 源质量、角向归一化、Up−In 符号 |
| spin-1/κ 的匹配方程 | `paper_all_component_matching.py` | m=1 的补充方程、角截断与相消精度 |
| h 与 κ0/κ2 径向积分及齐次补偿 | `paper_kappa_radial.py`、`lorenz_kappa.py`、`lorenz_chi.py` | 原 κ23 与本地 κ24/χ 的映射、角向交叉项、同一源跳跃 |
| 十个 tetrad 分量及投影 | `paper_full_tetrad.py`、投影诊断 | tetrad 缩放、rho/Sigma/Delta 因子、球谐而非椭球谐输出 |

Sam `.m:321–327` 明确使用负自旋的 Chandrasekhar 本征值。其 `.m:1355–1357` 的空间采样径向求解配置是 `NumericalIntegration`，MST/Monodromy 选项被注释；但更早求匹配幅度时有未显式指定 Method 的 TeukolskyRadial 调用。因此不能把整个原作者实现概括成“低频一律 MST”。

Sam `.m:218–221` 会把输入 a、r0 四舍五入到 3 位小数。要计算实际 a=0.8771530275949366，隔离驱动必须明确避免该配置读取精度损失，不能把 a=0.877 的结果误作相同参数。README 推荐至少 32 位精度，并指出高阶相消需要更高精度。

Conor 模块依赖 `OrbitalData`、`KerrOrbitalParameters`、`SWSHdecomp`、`SpinWeightedSpheroidalHarmonicsFT`、`NGrid`、`MSTMode`、`MST`，不是只安装标准 Teukolsky 包就能运行。其 RunFile 含个人绝对输出路径和 NotebookDirectory 设置；应先隔离驱动，不应直接执行整份文件。

## 两处不能忽略的分支差异

### 1. 四个副本会截掉 |m|=1 的 spin-1 偶极

在固定 Conor 提交中，定义是 `lmins0=Abs[mm]`、`lmins1=Max[1,Abs[mm]]`、`lmins2=Max[2,Abs[mm]]`。五个包含 `MMm1gridL` 的源码文件已全部核对：

| `.wl` 相对 Conor 根目录 | MMm1grid 的 If 行 | 条件 |
|---|---:|---|
| `GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl` | 2422 | `ll>=lmins2` |
| `GenerationCodes/Numerics-WorldTube/h1-radial-gen/MetricReconstructRadiative.wl` | 2431 | `ll>=lmins2` |
| `Ret-Numerics/BIGRUN-h1-radial-gen-session/MetricReconstructRadiative.wl` | 2431 | `ll>=lmins2` |
| `Ret-Numerics/h1-radial-gen/MetricReconstructRadiative.wl` | 2430 | **`ll>=lmins1`** |
| `old/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl` | 2422 | `ll>=lmins2` |

表中 4 个 `lmins2` 副本把 ell=1 的负自旋1径向数组置零；随后正自旋1由该数组构造，因此也为零。`Ret-Numerics/h1-radial-gen` 的条件正确。Sam 原 `.m:1209,1375` 均使用 `lmins1`，保留合法偶极。

这是真实的源码行为差异，**尚不是 2025 论文的错误证明，也不能仅凭目录名称认定某个 HDF5 就含有该错误**。当前分支同时存在正确和错误门槛。需要对 HDF5 做独立分量比对才能决定其可用范围。

### 2. 调试 loader 把两个 m=1 分量换成 m=2

`GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/LoadData.wl:42` 和 `old/` 对应副本包含：

```wolfram
miLoad = If[mi == 1 && MemberQ[{8, 9}, index], 2, mi];
```

这会在请求 m=1 的第 8/9 个分量时读取 m=2 的数据。该 loader 还有目录 a=.6 与文件参数 a=0 的不一致默认值。参考比较必须直接读取正确 HDF5 的 `/m_1` 数据集，不能通过这个 loader 生成真值。

## m=1 匹配方程：论文文字与公开实现不同

Sam `.m:884` 的基础系统由 l+l+ 与 m+m+ 的场值及导数条件组成。注释 `.m:890` 说 |m|=1 时需要 spin-1 分量补足；但**实际执行的 `.m:896` 是 ell=1 的 l−l− 值与导数条件**。混合分量版本在 `.m:893–895` 被注释。其余分量仍用于残差检查。

这些条件在完整一致的解上可以等价；有限角向截断与数值相消会影响条件选择的表现。因此它是值得逐条件对照的公开技术细节，不能仅依据注释认定作者代码使用混合分量，也不能未经传播就认为它解释百分数级通量偏差。

## κ 的正确映射与源匹配

作者的对角径向算子为

\[
L_\ell=\partial_r(\Delta\partial_r)+K^2/\Delta-\lambda_\ell,
\qquad L_\ell h_\ell=0,
\qquad L_\ell\kappa_{\ell\ell}=\tfrac12(r^2+a^2\Gamma_{\ell\ell})h_\ell.
\]

Sam `.m:1036` 分别解 `L κr0=h/2` 与 `L κr2=r²h/2`。`.m:1241–1263` 明确先由 trace 给 particular solution 的幅度，再加齐次 h 模式满足 κ 的值和导数跳跃。非对角项由 a²Γ 与本征值差得到，不能删去。

本地 mass-squared resolvent 的变量满足

\[
(\Box-\mu_s^2)h_{\mu_s^2}=16\pi T,
\qquad \kappa_{24}=-i\omega\partial_{\mu_s^2}h_{\mu_s^2}|_0,
\]

而原 2023 κ 的映射是

\[
\boxed{\kappa_{23}=\frac{i}{2\omega}(\kappa_{24}-\chi)
=\frac12\partial_{\mu_s^2}h_{\mu_s^2}|_0-\frac{i\chi}{2\omega}}.
\]

对应代码是 `report_paper_scalar_mapping.py:29–32` 和 `report_paper_kappa_radial.py:28`。χ 提供相关的齐次规范补偿。因此“作者逐模式构造 κ”与“本地质量变分”不能通过直接比较 κ24 数值判定相等或不同；必须使用上述映射、全角向交叉项、相同 retarded 边界与源匹配。质量变分在有限角截断时还会带入源投影的变分，完整角和才消去这些接触项。既往真空映射/对角径向检查不自动证明完整有源问题等价。

## 已实际执行原作者模块

本机 Windows Kernel 14.0 可运行，虽 stdout 为空，`Export` 输出成功。此次只原样提取：

- Sam `.m:43–59` 的微分算子定义，`.m:162–186` 的标量及 κ 方程。
- Conor 正确低阶版本 `.wl:64–88` 的 `ConstructSolution` 与 `ConstructSolutionhκ`。

驱动：`outputs/paper_metric_reference/radial_author_modules.wls`；结果：同目录 `radial_author_modules.json`。源码片段范围与 SHA-256 见 `code_audit.json` 和精简 manifest。

| 检查 | 实际结果 |
|---|---|
| 原作者标量算子与上式 L 的符号差 | 严格为 0 |
| 原作者 κ 算子与 `L κ−Σh/2` 的符号差 | 严格为 0 |
| 制造复数 trace 与 κ 跳跃，3 组 In/Up 基底缩放（含 1e±25） | 原作者 helper 的所有跳跃残差严格为 0 |
| particular κ 分别加任意齐次 h 后再匹配 | 最终 κ 严格不变 |
| 与本地 `PaperKappaRadial:105–110` 的 2×2 匹配代数比较 | 相对最大差 2.66e−16 |

这组检查排除基础径向方程和跳跃线性代数的约定不一致，**不是物理源幅度、完整度规或最终环境通量的独立验收**。没有运行作者全部角向匹配系统，也没有把制造数据检验称为 Kerr 原论文复现。

## 公开数据的读取边界

可用小文件包括 `GenerationCodes/Numerics-Asymptotics/h1-mmode-gen/h1mmodedat/data600/data600-80/data/h1_a0.6_rp8.0_l4_m{1,2,3}.h5`。文件名 l4 表示截断/输出包含至 ell=4，不能把它当作单一 ell=4 模式。

m=1 文件中 `/m_1/r_in` 为 101 点，r 从 8 向视界方向减少；`r_up` 为 3960 点，r 从 8 增至约 987.555。各 In/Up 分量分别是 `(5,101)` 与 `(5,3960)`，ell 行为 0..4。前九项是球谐分量；第十项 trace 必须单独确认：原 Sam `Calcs0components:2190` 返回 `hvec`，仍是椭球谐系数，需要乘 `Bmat0` 才得到球谐 trace。复数以 Re/Im compound dtype 存储。十个分量的 rho、Sigma、Delta 缩放必须按键名保留。未发现足以把数组绑定到上述某个 `.wl` 副本的生成提交元数据。

只为读 HDF5，在 `outputs/paper_metric_reference/python_deps` 隔离放置了 h5py 3.16.0，无依赖安装；未改动生产 `.venv`。未联系作者、未改动任何生产模块或原始下载文件。


## 最初 Sam 驱动的准备记录（后续运行状态见下节）

在 `outputs/paper_metric_reference/sam_L10_m1` 和 Windows `原作者代码对照_20260915/sam_L10_m1` 准备了完整物理参考驱动 `run_sam.wls`，由主任务启动。此报告撰写时仅做过 Kernel 解析检查，没有将“准备好”记录成已完成数值验收。

它原样保留 `.m:43–951` 的完整源匹配、`:954–1274` 的 trace/κ/自旋径向幅度，以及 `:1354–2248` 中物理替换与优化分量函数。删除的是原 CLI/目录读取、大网格和演示求值；配置和输出全部在隔离目录。没有改动下载原件。原 `.m` 匹配条件、相位、精度清理、ODE 和投影公式保持不变。

配置为 a=3/5、r0=8、m=1、L=10、prec=50、accgoal=24、nterms=12、kapord=10、inford=10、horord=8、rinf=4000、xhor=1/10000。只计算轨道两侧各一个与公开 HDF5 完全相同的 r：3.949919837349015、11.986125913940135。输出保留三类自旋扇区、十个原始数组，另单独输出球谐 trace，并在开始径向求解前先保存匹配跳跃/Bmat/本征值。

官方 BHPT 依赖由同轮独立目录的 `bhpt_paths.wl` 加载；没有全局安装。`01_matching.m`、`02_radial.m`、`03_fields.m`、主驱动分别成功解析 413、146、335、27 个表达式。源码边界和 SHA-256 清单为 `docs/environment_reproduction/sam_minimal_driver_manifest_20260915.json`。`status.json`、`run.log`、逐点 `sample_*.json` 用于区分实际运行阶段与完成结果。


## 原作者真实匹配数据的首次对照

Sam a=.6、r0=8、m=1、L=10 驱动已完成源匹配并进入径向阶段。新增 `src/report_author_sam_jump_comparison.py` 支持三种独立步骤：生成本地原 Sam 行选择的完整匹配、提取现生产 compact-current 的映射跳跃、读取已完成作者 JSON 做比较。它不轮询、不读取未完成的文件，不修改生产模块。

原 Sam `matching.json` 的数组只包含 ell=1..10，每行四项为 `[Pm1j0, Pm1j1, kappaj0, kappaj1]`；没有 ell=0 行。L=10 的受保护区只有 ell=1,2,3。现生产路径包含完整 current 的两种手性、κ23 的 χ 补偿与 angular source derivative，得到：

| ell | 4 项跳跃的行尺度相对最大差 | 本地 spin1 基底提取残差 |
|---:|---:|---:|
| 1 | 1.980e−14 | 5.631e−15 |
| 2 | 1.230e−12 | 9.564e−12 |
| 3 | 1.146e−13 | 6.133e−14 |

比较采用行中最大作者跳跃幅度作尺度，并另存绝对差；不是对对称性为零的 κ 项声称同样相对精度。这里的 spin1 基底提取是把实际规范矢量表示在固定的两个径向初值基底上，包含有限半径偏移后的二阶 Taylor 外推，残差已记录，不是对论文通量拟合参数。

这证明 **a=.6、r0=8 的实际原作者低阶源跳跃与当前生产映射相符**，比单纯内部闭合检查更强；它仍不证明后续径向积分、完整度规和环境通量相符。数据：`docs/environment_reproduction/sam_vs_production_jumps_ell1to3.json` 及相邻 `sam_local_production_jumps_ell1to3.json`。

此外，独立启动了实际偏差参数 a=0.8771530275949366、r0=20、m=1、L=10 的原 Sam 驱动，采样 r={2,10,30}。其隔离源码只有原第 219 行取消 3 位 Round；a 用精确有理数 8771530275949366/10^16 输入。启动 Kernel PID=1968，源码变更与配置见 `sam_actual_minimal_driver_manifest_20260915.json`。运行结果未完成前，不把这一步记作实际参数已验收。


### 实际偏差参数与独立行选择的完成结果

实际 a=0.8771530275949366、r0=20 的作者源匹配现已完成。其 ell=1,2,3 四项跳跃与生产映射的行尺度相对最大差依次为 **2.30e−14、3.36e−12、2.10e−13**，见 `sam_actual_vs_production_jumps_ell1to3.json`。这些是 source-jump 对照结果；实际参数的径向与 metric 仍单独验收。

a=.6、r0=8 的本地 L10/q32 完整独立匹配也已完成。隔离函数只改变 m=1 的补充方程选择为作者实际使用的 l−l− 值/导数条件，并保存同一组装矩阵重解原 mixed 控制。对作者真实匹配数据：

| ell | 原 Sam l−l− 行选择 | 同一矩阵 mixed 行控制 |
|---:|---:|---:|
| 1 | 2.33e−15 | 4.76e−13 |
| 2 | 4.39e−16 | 2.91e−14 |
| 3 | 4.01e−15 | 4.03e−15 |

这些是四项跳跃的行尺度相对差，说明补充方程选择不造成此配置的低阶明显偏差。低阶全分量残差的既定归一化最大值约 1.3e−11。没有隐藏角截断效应：ell=8,9,10 的跳跃差分别为 1.02e−8、1.70e−7、3.44e−6，整域未使用条件绝对最大残差为 2.81；它们不能与受保护 ell=1..3 的结果混为一谈。

完整记录为 `sam_vs_local_matching_L10_q32.json`、`sam_local_matching_L10_q32.json`。新增 harness 保存原函数与隔离改写的 SHA-256；生产 `paper_all_component_matching.py` 未修改。与作者代码源输入的一致性已获得真实原程序输出支持，百分数级环境通量偏差还需要继续由后续度规/源积分定位。


### 为定位径向阶段耗时新增的数值控制组

高阶驱动在径向符号级数构造阶段运行较久后，另启动 `sam_actual_L10_i6_h5`（PID 40264），保留原两组作对照。新组只把 inford=10→6、horord=8→5、rinf=4000→10000；方程、匹配条件、原物理函数和三个采样半径保持不变。这属于改变渐近边界的数值截断控制，不能仅称为缩短采样，也不能未经收敛比较就认定边界误差更小。

`01_matching.m`、`03_fields.m` 与先前实际参数组逐字节相同。`02_radial.m` 仅添加阶段和 kk/ell 日志；每一条原物理表达式保留。它会在匹配后 `DumpSave` 完整 Global 符号状态为 `matching_state.mx`，便于后续在相同背景、匹配参数下复用，省去再次匹配。来源行和改动哈希见 `sam_actual_i6_h5_driver_manifest_20260915.json`。新组三份源码已在同一个执行 Kernel 中以 Hold 方式完整解析通过；本文此处记录启动，尚非完成数值验收。


i6/h5 组进入分量采样后，新增精化组 `sam_actual_L10_i7_h6`（PID 37072），从已完成匹配的 `matching_state.mx` 恢复，不重算匹配。恢复前固定并校验 SHA-256，背景与配置 13 项检查全通过。仅覆盖 inford=7、horord=6、rinf=20000 和输出目录；带日志的 02/03 源码逐字节复用，采样点仍为 r={2,10,30}。其清单与恢复验证见 `sam_actual_i7_h6_driver_manifest_20260915.json`。精化组结果需完成后再与 i6/h5 定量比较，启动本身不构成收敛结论。


## 本轮最终完整度规执行结果

后续实际参数 a=.8771530275949366,rp20,m1,L10 的两组原Sam完整执行已成功：inford6/horord5/rinf10000 与 inford7/horord6/rinf20000。r2/10/30、保护球谐j1..3，十分量对本地最大相对L2差8.19e−7；两组作者边界精化最大相对变化3.34e−12。源码、配置、采样及严格范围见 [本轮总报告](AUTHOR_CODE_MODULE_COMPARISON_20260915.md) 和 [边界精化JSON](environment_reproduction/sam_actual_metric_boundary_convergence.json)。原inford10两组因符号级数成本终止；matching已完成且保留，但没有当作完整度规结果。下载原件和生产算法均未修改。
