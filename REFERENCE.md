# Research Resources

## 摘要描述

本文件是当前课题重要研究资源的渐进式入口，包括paper、个人note和已有代码。被收录表示值得优先查阅，不表示内容已经验证或仍然正确。先看索引；只有任务触发“深入查阅条件”时，才读取链接详情或原始资源。

除状态栏明确注明的例外（当前为 `P17`、`P18`、`P21`）外，paper均已在Zotero中定位。除状态栏明确记录已核对正文的条目外，概要只基于Zotero/arXiv元数据和abstract，不能用于支持具体公式、约定或适用范围；需要这些证据时，优先通过Zotero阅读正文。

## 具体内容

### Paper索引

| ID | 标题 | 支撑内容 | 何时深入查阅 | 状态 |
|---|---|---|---|---|
| `P01` | [Quasinormal-mode filters: a new approach to analyze the gravitational-wave ringdown of binary black-hole mergers](resource_details/papers.md#p01) | QNM滤波与次主导ringdown结构 | 使用滤波定义、稳定性或具体模态结论时 | Zotero已收录；filter相关正文已核对 |
| `P02` | [Black hole spectroscopy by mode cleaning](resource_details/papers.md#p02) | mode cleaning与黑洞谱学 | 使用Bayesian框架、likelihood或模型检验时 | Zotero已收录；metadata+abstract |
| `P03` | [Using rational filters to uncover the first ringdown overtone in GW150914](resource_details/papers.md#p03) | GW150914第一overtone | 使用Bayes factor、起始时间或振幅结果时 | Zotero已收录；metadata+abstract |
| `P04` | [The SXS Collaboration catalog of binary black hole simulations](resource_details/papers.md#p04) | SXS数值相对论波形目录 | 选择simulation、误差或参数覆盖时 | Zotero已收录；存在重复条目；metadata+abstract |
| `P05` | [A catalog of 174 binary black-hole simulations for gravitational-wave astronomy](resource_details/papers.md#p05) | 早期SXS数值相对论目录 | 使用目录参数、波形或precession结果时 | Zotero已收录；metadata+abstract |
| `P06` | [Probing Direct Waves in Black Hole Ringdowns](resource_details/papers.md#p06) | direct wave理论与可探测性 | 使用direct-wave定义、screening或高自旋行为时 | Zotero已收录；稳相相关正文已核对 |
| `P07` | [GW250114 reveals black hole horizon signatures](resource_details/papers.md#p07) | direct wave观测与horizon signatures | 使用GW250114数据、统计量或近视界解释时 | Zotero已收录；template相关正文已核对 |
| `P08` | [Analytic Black Hole Perturbation Approach to Gravitational Radiation](resource_details/papers.md#p08) | 黑洞微扰、PN与解析方法 | 使用Teukolsky、MST或PN公式和约定时 | Zotero正文相关方程已核对 |
| `P09` | [Gravitational radiation from plunging orbits: Perturbative study](resource_details/papers.md#p09) | plunge辐射与ringing激发 | 使用plunge轨道、频率或激发结论时 | Zotero正文的轨道与source附录已核对 |
| `P10` | [Ringdown of a postinnermost stable circular orbit of a rapidly spinning black hole](resource_details/papers.md#p10) | post-ISCO与高阶QNM激发 | 使用质量比、transition或高阶模结果时 | Zotero正文的transition与初值段落已核对 |
| `P11` | [The Transition from Inspiral to Plunge for a Compact Body](resource_details/papers.md#p11) | Ori-Thorne transition模型 | 使用transition scaling、方程或初始条件时 | Zotero已收录；transition相关正文已核对 |
| `P12` | [Gravitational Waves from a Compact Star in a Circular, Inspiral Orbit](resource_details/papers.md#p12) | 近ISCO圆轨道波形与flux | 使用表格、相对论修正或LISA估计时 | Zotero已收录；metadata+abstract |
| `P13` | [Numerical computation of second-order vacuum perturbations of Kerr black holes](resource_details/papers.md#p13) | Kerr Teukolsky双曲紧致化与数值演化 | 使用Eq. (22)--(25)、一阶化、RK4 constraint reset、特征边界或轴正则性时 | Zotero正文相关方程、Sec. V与Appendix C已核对 |
| `P14` | [Foundations of Direct Waves in Schwarzschild Ringdown](resource_details/papers.md#p14) | filtered DW的Green-function基础 | 判断anti-causal贡献、filter含义或BBH/Kerr外推边界时 | Zotero已收录；相关正文段落已核对 |
| `P15` | [The Direct Wave is Not a Meaningful Test of Horizon Properties](resource_details/papers.md#p15) | DW与horizon参数的NR检验 | 判断DW-horizon相关性、template系统误差或投稿新颖性时 | Zotero已收录；相关正文段落已核对；存在重复条目 |
| `P16` | [Hyperboloidal foliations and scri-fixing](resource_details/papers.md#p16) | hyperboloidal compactification与SCRI+边界 | 论证把null infinity纳入有限计算域、避免人工timelike outer boundary或解释scri-fixing时 | Zotero已收录；任务相关主张已核对 |
| `P17` | [On the Partial Difference Equations of Mathematical Physics](resource_details/papers.md#p17) | 双曲方程有限差分的CFL必要条件 | 从characteristic speeds和网格尺度设计显式时间步或讨论收敛必要条件时 | 英译原文与书目信息已核对；未逐式审阅全文 |
| `P18` | [Spin-s Spherical Harmonics and ð](resource_details/papers.md#p18) | spin-weighted harmonics与edth算符 | 推导固定m场在轴上的正则幂次、angular operator cancellation或SCRI+模态投影时 | 原文与书目信息已核对；未加入Zotero |
| `P19` | [Notes on the integration of numerical relativity waveforms](resource_details/papers.md#p19) | 从 $\psi_4$ 恢复 strain、fixed-frequency integration | 设计 $\psi_4\to h$ 二次积分、低频 cutoff、time window 或 sampling 检查时 | arXiv v3 正文相关公式与数值讨论已核对 |
| `P20` | [How much of the black hole ringdown is sourced within the light ring?](resource_details/papers.md#p20) | Kerr粒子源的LR内外切分、filter和可测源区 | 判断NH公式与空间源区是否等价、mask及时间补偿的含义时 | Zotero正文及arXiv v1相关内容已核对；未复现其数值 |
| `P21` | [Complex frequency evolution of direct waves from binary black hole mergers](resource_details/papers.md#p21) | NR/CCE动态频率、Kerr pole–zero结构与校准 | 比较晚时极限、News/strain、filter选择和horizon推断时 | Zotero按ID/标题未定位；arXiv v1相关正文已核对 |

详细元数据和abstract概要见 [paper条目](resource_details/papers.md)。

### 本地与历史资源索引

| ID | 标题 | 类型 | 用途 | 访问与可信状态 |
|---|---|---|---|---|
| `C01` | [Direct_Wave_in_Higher_order/Codes](resource_details/local_resources.md#c01) | code collection | 历史代码入口 | T009 另授权主 agent 阅读 errorAnalyse、mma-expansion 框架；coding 子 agent 禁读；历史内容不自动可信 |
| `H01` | [Paper/main.tex](resource_details/local_resources.md#h01) | manuscript draft | 历史文章记录 | T004 已审计数值演化附录；发现约定缺口；未批准 |
| `H02` | [Paper/ExciseOfDW.md](resource_details/local_resources.md#h02) | note | 历史想法与推导记录 | 未读取；不保证正确 |
| `H03` | [direct_wave_quadratic.tex / PDF](resource_details/local_resources.md#h03) | manuscript draft | 非线性direct-wave历史记录与广义Lu模板 | T010已核对owner提供PDF第4页p模板；未验证全文 |
| `H04` | [Groupmenting824.md](resource_details/local_resources.md#h04) | meeting note | 组会历史记录 | 未读取；不保证正确 |
| `N01` | [Mathematica MCP usage guide](resource_details/local_resources.md#n01) | note | MMA调用与故障处理 | 已核对全文；环境信息需复查 |

本地资源详情只记录路径、用途和访问边界，见 [本地与历史资源](resource_details/local_resources.md)。历史材料不能直接支持可信物理结论；真正可信的内容依赖后续正文核对、明确推导或实际执行检查。

## 总结及下一步预告

当前已登记21篇paper、1个受限代码集合、4份历史文档和1份MMA指南。后续任务先从本索引定位资源，再按条件加载详情；不要一次性读取全部条目。
