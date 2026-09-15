# 官方 Mathematica BHPT 依赖及函数对照（2026-09-15）

已直接运行官方 Mathematica BHPT 函数，得到与当前 `pybhpt 1.0.0 / AUTO` 独立吻合的轨道、角向及径向参考。所有 11 组、66 个径向采样位置上的 `R` 和 `R'` 均有限；保持原始单位传输归一化，未拟合相位或振幅。

## 依赖来源与加载

本机 `/home/ljq` 有界文件搜索未找到可直接复用的完整官方 Mathematica Teukolsky/SWSH/KerrGeodesics 包；Conor 仓库的 `SpinWeightedSpheroidalHarmonicsFT.m` 是作者自定义模块，没有冒充官方 paclet。现有 Julia GSN 项目不提供这套 Mathematica 加载接口。

从官方 Black Hole Perturbation Toolkit GitHub 下载了三个小型源包，压缩后总计 848204 字节，均保留 MIT LICENSE。

| 官方包 | paclet 版本 | 固定 commit |
|---|---|---|
| [Teukolsky](https://github.com/BlackHolePerturbationToolkit/Teukolsky/tree/a7a508b0799ca8ac48a42fa82a9223df962499b0) | 1.1.1 | `a7a508b0799ca8ac48a42fa82a9223df962499b0` |
| [SpinWeightedSpheroidalHarmonics](https://github.com/BlackHolePerturbationToolkit/SpinWeightedSpheroidalHarmonics/tree/49f9b87d786f1b5323678fa83dce8b7e4aeea0ec) | 1.1.0 | `49f9b87d786f1b5323678fa83dce8b7e4aeea0ec` |
| [KerrGeodesics](https://github.com/BlackHolePerturbationToolkit/KerrGeodesics/tree/c7c57a936236ddcce98963481440bca0c433b3c3) | 0.9.0 | `c7c57a936236ddcce98963481440bca0c433b3c3` |

这些是本次检索的官方固定版本。提交日期分别为 2025-06-26、2026-09-14、2026-09-15，**不是已确认的原论文 2025 年历史依赖锁定文件**。

`Teukolsky/Kernel/Teukolsky.m` 检查 SWSH>=1.0.0、KerrGeodesics>=0.9.0；MST 位于该包自身的 `Kernel/MST`，不需另装 `MST` 包。使用 Windows 原生 Mathematica 14.0.0，仅通过 `PacletDirectoryLoad` 临时加载 WSL 中的 `outputs/paper_metric_reference/wolfram_dependencies/`，没有 `PacletInstall`、全局安装、虚拟环境变更或生产代码修改。

已为 Sam 原 `.m` 驱动提供共享加载文件 `原作者代码对照_20260915/bhpt_paths.wl`。本环境需要转义反斜杠 UNC 路径；以正斜杠写 `//wsl.localhost/...` 的初次探测未正确加载，因此最终驱动使用前述 loader，主 `.wls` 放在 Windows 本地工作目录。

## 官方函数和参数

直接调用 `KerrGeoOrbit`、`KerrGeoFrequencies`、`SpinWeightedSpheroidalEigenvalue`、`SpinWeightedSpheroidalHarmonicS` 及其导数、`TeukolskyRadial` 及其导数，没有重写这些函数的物理公式。

- Sam 参考参数：`a=3/5, rp=8, ell=2, m=1, s=-2,-1,0,1,2`；半径为 `3,8,20`。
- 实际偏差诊断参数：`a=8771530275949366/10^16, rp=20, ell=1, s=-1,0,1, m=+1,-1`；半径为 `3,20,40`。
- 每组 `omega=m Omega_phi`，频率从官方轨道包获得。角向在 `theta=pi/4,pi/2,3pi/4, phi=0` 采样，`SphericalExpansion` 使用 30 项。
- 径向明确使用官方 `Method -> "MST"`、`WorkingPrecision -> 40`、`AccuracyGoal -> 25`、`PrecisionGoal -> 25`。每个 case 的时间上限为 120 秒，11 组均完成。
- 同时导出 In/Up 的本征值、重整化角动量 `nu`、入射/透射/反射振幅和逐点 `R,R'`。`Transmission=1`。

每个 JSON 包含便于 Python 读取的数值，以及保留 Wolfram 任意精度数字的 `raw_input_form`。这里的 40 位工作精度是计算设置，不是对整个物理流水线误差的承诺。

实际 `rp=20` 官方轨道结果为：

```text
Omega_phi = 0.0110717605820720444302713976564...
E         = 0.975631952782719837165183191498...
Lz        = 4.71086066364514414560491798055...
u^t       = 1.08286701398582140793851010226...
```

实际低自旋径向的官方 MST 结果：

```text
nu(s=0,  ell=1,m=±1) = 0.9996911721643061...
nu(s=±1, ell=1,m=±1) = 0.9996175886816968...
```

这与此前隔离修正的低频 `nu` 及 MST 锚点验证相符；这里使用的是独立官方 Mathematica 实现。

## 与现有 AUTO 的数值比较

[逐项 JSON](wolfram/official_bhpt_AUTO_comparison.json) 保存全部比较；[比较脚本](wolfram/compare_official_bhpt.py) 只读取官方输出，并调用现有 `pybhpt`，不改运行环境。

| 量 | 11 组最大差异 |
|---|---:|
| 原始 `R` 相对差 | `3.0065e-14` |
| 原始 `R'` 相对差 | `3.3356e-14` |
| 角向本征值绝对差 | `1.3172e-13` |
| `S(theta,0)` 绝对差 | `5.8009e-15` |
| `dS/dtheta` 绝对差 | `1.3989e-14` |
| 官方三个半径的加权 Wronskian 相对漂移 | `4.4710e-16` |
| 缩放不变 Green 函数样本相对差 | `2.1763e-14` |

加权 Wronskian 使用 `Delta^(s+1) (R_In R_Up' - R_Up R_In')`。Green 样本取 `R_In(r_min) R_Up(r_max) / W(r_mid)`。这是物理解及归一化相容性的附加诊断；原始函数比较已经无需相位或振幅调整而一致。角向相对归一化附加检查避开零节点，主要误差表使用绝对误差。

这组独立参考不支持“这些实际低自旋径向函数、角向本征值或角向导数本身产生约 32% 差异”的解释。范围限于这些自旋、频率、角度和半径；不是对未采样的近视界端点、远域、较高 ell 或最终环境源/通量的完整证明。后续应继续把作者度规、补全和环境源逐层对照。

## 复现与文件

- [下载及校验脚本](wolfram/fetch_official_dependencies.py) 按固定 commit 和 tar SHA 获取这三个依赖到 ignored outputs；已有文件必须相同，拒绝覆盖不同内容。
- [官方函数驱动](wolfram/bhpt_official_probe.wls) 与 [临时加载器](wolfram/bhpt_paths.wl) 供 Windows Mathematica 原生运行。加载器目前含本机已验证的 UNC 根路径，换机器时修改这个路径。
- [单独轨道导出脚本](wolfram/bhpt_orbit_export.wls) 修正首次导出将 `FourVelocity` 的列表误当 association 的格式问题；计算函数及径向输出未改变，最终轨道 JSON 完整。
- [manifest](wolfram/official_bhpt_manifest.json) 记录依赖 URL/commit/源包哈希和所有结果、脚本哈希。
- 完整官方原始函数结果位于 `outputs/paper_metric_reference/wolfram_dependencies/official_results/`，三个依赖自身另有逐文件源哈希清单。

在 repo 根目录运行：

```bash
.venv/bin/python docs/environment_reproduction/wolfram/fetch_official_dependencies.py
.venv/bin/python docs/environment_reproduction/wolfram/compare_official_bhpt.py
```

官方 `.wls` 已在本机以 `D:\mathematica\mma\WolframKernel.exe -noprompt -script <本地脚本路径>` 成功执行。标准输出为空，因此以导出的 JSON 和进度文件核验结果。

版本控制中的三份官方驱动/loader副本仅统一LF行尾和末尾空行；实际执行过的原字节另存`outputs/paper_metric_reference/wolfram_dependencies/executed_script_bytes/`，两种SHA与非物理变更说明见manifest的`archival_script_normalization`。
