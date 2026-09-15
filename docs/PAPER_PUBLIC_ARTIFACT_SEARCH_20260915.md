# Dyson et al. 原论文公开产物检索（2026-09-15）

这次找到了**第一作者及 Lorenz 重构作者公开的度规代码仓库**，并保存了原论文的 TeX、独立图 PDF 和作者接受稿。尚未定位到能够确认为这篇 2025 论文完整环境标量求解器的公开发行版，也未定位到其通量曲线的原始数值表。后两项是有范围的检索结果，不表示文件不存在。

机器可读来源、版本、请求结果和全部 49 个保存文件的 SHA256 在 [manifest](PAPER_PUBLIC_ARTIFACT_SEARCH_20260915_manifest.json)。材料保存在 `outputs/paper_original_reference/`，该目录为忽略输出，不修改生产代码，不改虚拟环境。此次未联系作者，也未下载受限附件。作者 Lorenz 仓库的代码下载和逐模块执行审查由另一子任务负责，本报告仅锁定其来源。

## 1. 原论文版本与实际保存的材料

- [arXiv:2501.09806](https://arxiv.org/abs/2501.09806) 当前只列出 v1，提交于 2025-01-16 19:32:12 UTC。
- [固定 v1 源码包](https://arxiv.org/src/2501.09806v1) 已再次下载验证，与原本地副本逐字节相同。这里的“源码包”是论文排版源文件，并非数值计算源代码。
- [White Rose 作者接受稿记录](https://eprints.whiterose.ac.uk/id/eprint/227286/) 挂有 [接受稿 PDF](https://eprints.whiterose.ac.uk/id/eprint/227286/8/Self_force_and_environments%20%282%29.pdf)，文稿日期为 2025-06-03。PDF 包括机构封面及正文、补充推导，共 14 页。机构标记作者接受稿、CC BY 4.0。新下载也与本地副本逐字节一致。
- [MPG 机构 API](https://pure.mpg.de/rest/items/item_3633555_3) 记录更新于 2025-06-05，提供 [公开预印本 PDF](https://pure.mpg.de/rest/items/item_3633555_3/component/file_3633556/content)。下载后 MD5 与机构元数据一致。

| 文件 | 字节数 | SHA256 |
|---|---:|---|
| `2501.09806v1.tar` | 3381880 | `c07028320ab015cb3dc8c113614a9ad785de77169e421c608a9f72efd0039665` |
| `accepted_2025.pdf` | 2537467 | `69a6ad0fc933bae8774b1673e001016d2a4c82c8c61984665d118c368158099c` |
| `arxiv_v1_MPG_preprint.pdf` | 4110090 | `6f41f28f1e89e518b4bd252390d017e552d80bbd762586aa7c2d64a141eb6486` |

`source_download_verification_20260915.json` 保存前两项重新下载的状态和哈希核验。

安全展开源 tar 后，其文件清单为：

```text
main_PRL.tex
main_PRL.bbl
ref.bib
references.bib
Figures/BosonEMRIsFieldPlots.pdf
Figures/Flux_Inf_Hor_TwoPanels_2.pdf
Figures/Flux_l_mode_convergence_inf.pdf
Figures/Flux_l_mode_convergence_particle.pdf
Figures/RelDiffFlux_Inf_Hor_TwoPanels.pdf
Figures/ResonanceBosonEMRIsFieldPlot.pdf
Figures/plotTotal_onepanelSM_Ratio.pdf
```

这个包中没有 `.m/.wl/.nb/.py` 等数值求解程序，也没有 CSV、HDF5 等原始通量表。独立图 PDF 是作者原图，可用于数字化和校验图线；从矢量路径或图像提取的点仍然是图像数字化结果，不能称为作者原始数值数据。完整清单与每张图的哈希已保存。

对源 TeX、接受稿可提取文字检索 `github`、`zenodo`、`figshare`、data/code availability、available upon/on request、publicly available，未发现环境程序或数据的仓库地址。不能据此推断正式期刊版是否另有数据声明。

## 2. 期刊补充材料的可访问范围

[APS 正式页面](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.134.211403) 给出 PRL 134, 211403，2025-05-28 发表，2025-04-25 接受；当前页面将补充材料标为 **Supplemental Material (Subscription Required)**。补充入口未提供本次能够验证并下载的独立文件，正式正文同样显示授权要求。

[哥本哈根大学记录](https://researchprofiles.ku.dk/da/publications/environmental-effects-in-extreme-mass-ratio-inspirals-perturbatio/) 列出提交稿和出版版本；其 [出版 PDF 地址](https://researchprofiles.ku.dk/files/454207414/PhysRevLett.134.211403.pdf) 在本次抓取返回 403。

MPG API 更明确地区分了附件：`file_3633556` 预印本为 `PUBLIC`，`file_3653088` 出版版为 `AUDIENCE`。没有尝试读取受限附件；元数据中仅列这两项，没有单独数据/程序附件。因此：本次可以审查原 arXiv/接受稿所包含的补充推导，**不能声称已完整审查 APS 单独补充文件或正式版数据声明**。

## 3. 作者公开 Lorenz 度规代码：实质性新来源

| 仓库 | 固定 commit | 来源关系与限度 |
|---|---|---|
| [ConorDyson/KerrLorenzMSF](https://github.com/ConorDyson/KerrLorenzMSF/tree/6b469f0835e42d3fa531b4997d7a912257119f60) | `6b469f0835e42d3fa531b4997d7a912257119f60` | 第一作者公开 fork，建于 2026-01-15，父仓库为下项；晚于原论文，不能自动认定是论文当年的完整运行版本。 |
| [srd24/KerrLorenzCirc](https://github.com/srd24/KerrLorenzCirc/tree/5c1b42793ff893fd0c65a5b48ceecc02d34c5e35) | `5c1b42793ff893fd0c65a5b48ceecc02d34c5e35` | Sam Dolan 的公开 Kerr 圆轨道 Lorenz 度规重构代码，仓库创建于 2025-04-07。提供可信的作者代码独立对照。 |

两个仓库的 repo/commit/recursive-tree API 原文已保存。前者 tree 有 2183 项，后者 66 项。代码下载、许可证处理、符号及数值接口审查另见主任务的 Lorenz 作者代码审查产物。

这改变了“只能根据论文自行重写度规重构”的局面：可以直接对照作者公开实现。但是，**作者提供的 Lorenz 重构依赖，不等于这篇论文完整的云、协变源、投影、质量标量 Green 函数和最终通量流水线**。对 2025 论文具体补全约定、归一化和环境程序关联仍需逐项核对。

## 4. 其他作者代码和数据索引的检查

2026-09-15 的公开仓库单页列表均不足 100 项，已完整保存：ConorDyson 4 项、thomasspieksma 7 项、richbrito 5 项、srd24 4 项。这些是公开列表，不包括私有仓库。AEI 的 [Maarten van de Meent 公开 GitLab](https://git.aei.mpg.de/mmeent) 经 API 返回一项课程材料项目。

以下相近项目容易被误当为本论文代码，因此保存 README 并锁定 commit：

| 项目 | commit | README 对应的实际研究 |
|---|---|---|
| [thomasspieksma/GrAB](https://github.com/thomasspieksma/GrAB/tree/648ec51659e327f49125f3897566bd9c6cd9282d) | `648ec51659e327f49125f3897566bd9c6cd9282d` | 引用 2305.15460、2403.03147、2407.12908；引力原子双体共振/电离程序及预计算数据。没有声明为 2501.09806 的完整相对论环境求解器。 |
| [richbrito/gw_superradiance](https://github.com/richbrito/gw_superradiance/tree/6f482a2adecd7dad0a36546397e116e89915ada8) | `6f482a2adecd7dad0a36546397e116e89915ada8` | 标量云引力波发射，主要用于 1706.06311、1706.05097；不是这里的轨道诱导标量通量表。 |
| [richbrito/Tidal_Grav_Atoms](https://github.com/richbrito/Tidal_Grav_Atoms/tree/b96b71348dffd61253fa482e7cf3586a2e4ddab1) | `b96b71348dffd61253fa482e7cf3586a2e4ddab1` | 2410.00968 引力原子潮汐 Love 数。 |
| [ConorDyson/SpiralDensityWavesInKerr](https://github.com/ConorDyson/SpiralDensityWavesInKerr/tree/78de2c1ae06abdc563c5995504e53a62908e2bf0) | `78de2c1ae06abdc563c5995504e53a62908e2bf0` | Dyson 与 D'Orazio 另篇螺旋密度波论文的推导 notebook。 |

这些仓库可能提供相邻模块的参考，但不能使用其中数据冒充本论文原始通量。

索引检索结果：

- Zenodo：以论文 DOI、arXiv 编号、完整标题短语分别检索，均 0 命中。
- DataCite：以 DOI 和 arXiv 编号分别检索，均 0 命中。
- Figshare：按 `resource_doi=10.1103/PhysRevLett.134.211403` 查询，返回空列表。保存文件名 `figshare_title_search_20260915.json` 沿用最初命名，**实际查询是 DOI 过滤**，准确 URL 在 manifest。
- Crossref DOI 元数据 `relation={}`，没有关联数据 DOI。
- GitHub repository metadata 搜索 `2501.09806` 为 0；`environment Kerr` 的 5 项为无关项目。没有进行需要登录的全站代码搜索。
- 网页检索论文 DOI/arXiv、原通量 PDF 文件名、作者与代码词的组合，没有新增可确认的原始数据来源。
- 已有 Spieksma 2025 博士论文可提取全文/相关章节检索没有找到环境程序公开地址；这属于补充交叉检查，不是原始数值数据。

所有“0 命中”只适用于这些索引和查询；文件可能以其他标题存档、未被索引或未公开。

## 5. 对当前逐模块对照的实际含义

| 模块/对照层次 | 本次取得的原始来源 | 可作出的比较 |
|---|---|---|
| Kerr Lorenz 度规及低自旋重构 | 作者公开固定 commit 代码 | 可进行真正独立的作者实现对照，优先核查四维分量、基、源幅值及补全项。 |
| 背景云、二阶协变导数与环境源 | 原论文/接受稿补充推导 | 尚无确认的完整原环境程序；依然要标为公式复现与独立算法测试。 |
| 环境质量标量径向/投影/积分 | 原论文数值流程说明 | 尚无原论文运行配置和完整数值程序，不能把本地实现之间的一致性称为与作者代码一致。 |
| 最终无穷远/视界通量 | 原 TeX 附带独立图 PDF | 可做带数字化不确定度的图曲线比较；尚不能做作者原始数值表逐点比较。 |

因此本次最可执行的新增路线是用锁定的 Lorenz 作者代码独立核对度规与源的输入，然后继续排查环境部分差异。原始通量表和完整环境求解器的公开位置仍未确定。
