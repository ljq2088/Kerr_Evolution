# Zotero“双曲框架”来源与阅读记录

读取日期：2026-09-10。收藏夹键：`XVNLY9AK`，名称“双曲框架”。
已读取该收藏夹全部条目元数据；原始条目标题、条目键、citation key、DOI、链接的简表保存在
`reference/zotero_collection.json`。Zotero 原库未被修改。

## 直接用于本次实现

| 文献 | Zotero 条目 / 附件 | 阅读与用途 |
|---|---|---|
| Macedo & Zenginoglu, Hyperboloidal Approach to Quasinormal Modes (2025) | `8AW4IT3C` / `SBPNQGZP` | 读取索引正文；第 3.2 节的高度函数、场正则化、紧致化与边界思想。它是框架入口，不能单独当作完整 Kerr 方程。 |
| Ripley, Computing the quasinormal modes and eigenfunctions for the Teukolsky equation using horizon penetrating, hyperboloidally compactified coordinates (2022) | `R7UUGT4R` / `RCDYXJNN` | 读取索引正文，重点核对第 2 节式 (5)–(11)。本项目与式 (11) 的场自旋权重零情形一致。 |
| Zhen-Tao He, Hyperboloidal framework in black hole dynamics | `KF65HIFK` / `ZYBLE3FY` | 读取讲义索引正文，重点第 11–14 页关于高度函数、Kerr 延伸的讨论；用 Ripley 原方程核对，而不直接复制讲义抽取文本中的度规分量。 |

综述引用的 Macedo (2020), *Hyperboloidal framework for the Kerr spacetime*,
[arXiv:1910.13452v2](https://arxiv.org/html/1910.13452v2)，也已核对第 4.1 节及附录 C。
本次收藏夹清单未显示它是其中的独立条目；该文通过公开原文补充阅读，未导入 Zotero。
本项目的 `sigma=r_plus/r` 属于 radial function fixing minimal gauge。

## 收藏夹内其他文献的处理

- Ripley 等 (2021)，*Numerical computation of second order vacuum perturbations of Kerr black holes*，`6TFRLV7N`：已查看元数据与摘要。附件 `AT863U9B` 本地全文索引接口返回 404；本次未将其作为系数核对依据。
- Zenginoglu (2008)，*Hyperboloidal foliations and scri-fixing*，`THVCDNRH`：已查看元数据与摘要。附件 `X9JMDTQQ` 的全文索引接口返回 404；实现使用上表中已读取的 Kerr 原方程。
- Leaver (1985)，*An analytic representation for the quasi-normal modes of Kerr black holes*，`VYSELPRX`：元数据与摘要检查，用于识别频域方法背景；本次未实现 Leaver 连分式。
- *Binary black hole coalescence in the large-mass-ratio limit: the hyperboloidal layer method and waveforms at null infinity*，`43PJUEFC`：元数据检查，属于相关波形工作。
- *Gravitational radiation from Kerr black holes using the Sasaki-Nakamura formalism: waveforms and fluxes at infinity*，`62UDB62L`：元数据与摘要检查，采用不同扰动变量，本次未调用该 formalism。
- *Bondi-Sachs Formalism*，`ALTDVDYY`，及 *Bridging time across null horizons*，`HLY4BCNT`：元数据检查，未用于具体系数。
- `H28HGKBA` 的元数据标题是 *Identify Web-page Content meaning using Knowledge based System for Dual Meaning Words*，与任务不符；没有据此引用物理结论，也没有改动该条目。

`BK7IBKZU` 是同收藏夹内另一个 He 讲义独立附件，本次读取的是关联论文条目的 `ZYBLE3FY`。
条目键（如 `R7UUGT4R`）与 BibTeX citation key（如 `ripleyComputingQuasinormalModes2022`）不同。
Zotero 的 BibTeX 导出接口两次返回 HTTP 500，因此 `reference/references.bib` 根据已读取的条目元数据生成，
保留现有 citation key；不将其描述为 Zotero 原生成功导出。
