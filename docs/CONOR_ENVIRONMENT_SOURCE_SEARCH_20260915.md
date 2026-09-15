# Conor 仓库环境程序全文检索（2026-09-15）

**在所下载固定版本的全部可读文本中，未找到可识别为原论文完整云背景、有质量 Klein–Gordon 环境源或标量通量求解器的程序。找到了真实的 `boson` 度规数据输出痕迹，不能将其忽略，但它们并不构成完整环境程序。**

范围是 [ConorDyson/KerrLorenzMSF 固定 commit](https://github.com/ConorDyson/KerrLorenzMSF/tree/6b469f0835e42d3fa531b4997d7a912257119f60) 的本地 2111 个文件。其中 2069 个文本文件共 47,347,857 字节全部搜索，包括 `.m/.wl/.wls/.nb/.md/.txt/.log` 及可读的无扩展名文件；排除 14 个 HDF5、1 个 PDF、27 个 `.DS_Store` 二进制文件。每个搜索文件均计算 SHA256，并以 Git blob SHA 验证与该 commit 的 GitHub tree 相同，没有仅搜索根 README。一份 `.log` 含 NUL 填充，已去掉填充后搜索其可读文本，并在 manifest 单独注明。

[manifest](CONOR_ENVIRONMENT_SOURCE_SEARCH_20260915_manifest.json) 包含完整文件清单、哈希、排除原因、全部模式、34 个候选文件的 362 处上下文，以及人工分类。

| 潜在线索 | 核对后的含义 |
|---|---|
| `old/metric_reconstruction_calc_radiative.nb:35807` 及副本：`lm_in_-boson-60-060-001.bin` | 已保存的输出单元格，路径指向 Sam 的 KerrLorenz **度规模态**数据；附近还列出 `lm_in/up_s2/s1/s0`。是历史 boson 标签，未见其后有环境云求解入口。 |
| `old/metric_reconstruction_calc_static.nb:10633` 及副本：`jumps_completion_-boson-00-060-000.mx` | 静态度规跳跃/补全输出，与环境程序可能使用度规数据相容，但不能据此认定是环境标量数据。 |
| `CalcStaticRadialFunctions.nb`、`CalcStaticTensorMode.nb`、`Jigsaw.nb` 中的 scalar / scalar sources | 上下文为 free scalar、trace、恢复 Lorenz gauge 和 completion，属于度规重构的自旋零辅助部分。 |
| `Derivation/MRR_dev.nb:333` 的 “Spin 0 solver” | 进一步检查 `teukP0/teukh0`：无质量自旋零 Teukolsky、迹和 κ 方程；不是含云质量参数的环境 KG 方程。 |
| `old/metric_reconstruction_calc_radiative.m:474` 起的 Flux tests | 直接调用 `TeukolskyPointParticleMode[-2,...]`，比较的是自旋二引力通量。`flux.cpp` 的零星引用位于旧度规验证 notebook。 |
| `alpha` | 已检查的活跃包中是线性匹配振幅，旧式表达式还定义 `AlphaSq=a^2-a*m/omega`，不能当作云耦合 `M mu`。 |
| `Ret-Numerics/h1-mmode-gen/Evaluation$h1Ret.nb:113` 等 | 路径明确含 `Self-Force/Gravitational-SF`、`SecondOrderSource`、`h1Ret/h1Punc`。仓库名 MSF 在已检索文本中未正式展开，不能凭缩写将其解释为 massive scalar。 |

检索计数：精确 `2501.09806` / 正式 DOI 0 处，massive / Klein / Gordon / Leaver / hydrogenic / quasibound 0 处，cloud/ScalarField/environment 等特定词 0 处；`boson` 36 处，均为上述两个旧 notebook 及重复副本中的保存路径。泛 `scalar` 250 处、`flux` 68 处、自力相关路径 8 处，已按上下文分类。面向常见文本及 Mathematica `SuperscriptBox` 的 μ² 模式无命中；这只是额外线索，不能单独证明某个方程不存在。

原先无边界匹配 `211403` 得到的五处候选，是长浮点常数和计时小数中的偶然子串；改用精确论文标识后均排除。补查新增的 1042 个 `.log` 和可读无扩展名文件没有增加相关候选。42 个 notebook 按哈希只有 12 种不同内容，另检查其章节标题，未发现隐藏的独立云/环境求解章节。

结论限于这个固定版本已下载文件的完整可读文本。没有执行 notebook，没有解码 `CompressedData`，也没有读取 HDF/二进制或作者本机被引用但未公开的文件。因此应继续表述为“已找到作者度规实现及其 boson 数据标签，完整原论文环境标量流水线仍未定位”，而不是“作者没有这套代码”。此次没有启动数值计算或修改任何作者源文件。
