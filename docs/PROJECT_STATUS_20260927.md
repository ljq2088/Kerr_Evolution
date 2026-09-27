# 项目状态与复现入口（2026-09-27）

本分支 `codex/reproduce-kerr-environment` 保留原来的双曲切片演化，并继续实现 Kerr 标量云的环境扰动、Lorenz 度规重构、通量和场图比较。后者是固定圆轨道的频域受迫问题，不是把原来的无源标量脉冲演化改名。

## 已完成结果

参数 M=1、a=.88、mu M=.3、(ell_c,m_c,n_c)=(1,1,0)、r0=20M。完整 ell<=5 的18个非静态有效轨道通量通道：

| 单位 epsilon²(Mc/M) 的有效能量通量 | 本地 | Li 原矢量图读数 | 相对差 |
|---|---:|---:|---:|
| 无穷远 | 1.8978444432481997e-5 | 1.8943078476349247e-5 | +0.186696% |
| 视界 | -6.828254573842368e-5 | -6.615228308272457e-5 | +3.220241% |

[完整计算记录](li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512/execution.json)保存逐模结果、输入来源和比较口径。这里的参考来自图中矢量曲线，不是作者原始数组；同参数有限分辨率结果接近参考，并不证明所有截断已收敛。论文归一化 epsilon² zeta² 与上表相差 alpha^6，其中 zeta=alpha^3 sqrt(Mc/M)。

旧 Dyson 原视界曲线后来被修正，不能继续把旧曲线的约32.7%偏差当作本程序的未解释误差。[参考修正审查](KERR_HORIZON_REFERENCE_CORRECTION_20260916.md)。

## FIG.1：仍在后台计算

目标是 Dyson 原文 r_p=3.5M 的赤道面及子午面场图。采用精确 a=.88 的云、度规及传播算子；正负 m、静态驱动和束缚模式均保留，ell=2..12 共88模。

本次发布的是70个已完成模态的冻结快照，未发布完整 FIG.1，也不把未完成模态补成零。以 [快照清单](figure1_alignment_20260926/snapshot_20260927.json) 中的 UTC 时间、完成清单和 SHA256 为准。远端的快照不会随本机任务自动更新。

[后台方案与限制](figure1_alignment_20260926/README.md)。最终程序会相干合成复场后取绝对值，比较 ell<=5、8、12，并使用原图固定色标，无拟合幅度或相位。

场图采用固定圆轨道。辐射通量没有反馈到 r_p(t)；这与论文生成场分布图的设置相符。Li 的长期轨道相位偏移属于另一个能量平衡计算阶段，当前场图任务未实现该阶段。

## 模块对应与关键约定

| 模块 | 实际职责与易混淆处 |
|---|---|
| [li_leaver_cloud.py](../src/li_leaver_cloud.py) | 求复频准束缚谱及150/300项径向级数；本身不做质量归一化。 |
| [li_normalized_cloud.py](../src/li_normalized_cloud.py) | 保留复频空间剖面，显式冻结时间增长；有限 BL t=0 Killing 质量归一化。不是把 KS 质量换个名称。 |
| [li_separated_source.py](../src/li_separated_source.py) | 由加权球谐 tetrad 系数组装半解析源；p<=12、64位 Fourier 系数，最终收缩仍为双精度。 |
| [li_source_terms.json](../src/li_source_terms.json) | 符号生成的86个印刷展开单项式及3个由协变恒等式要求补回的单项式；保留印刷式控制开关。 |
| [environment_lorenz_mode.py](../src/environment_lorenz_mode.py) | 本地独立 Lorenz 重构；椭球输入 L<=20 与球谐输出 j<=18 是两种截断。 |
| [environment_static_lorenz.py](../src/environment_static_lorenz.py) | 静态匹配和 completion；保留边界/跳跃检验限制。 |
| [li_order4_green.py](../src/li_order4_green.py) | 实正频率的四阶 In/Up 边界与径向 Green 函数。 |
| [li_field_green.py](../src/li_field_green.py) | 用 (omega,m)→(-omega,-m) 共轭关系处理负频率延迟解，不错误地选择入射支。 |
| [environment_response.py](../src/environment_response.py) | 在原源积分分段内重构连续场，从两端分别累计 Green 积分，减少相减消失。 |
| [report_li_aligned_flux.py](../src/report_li_aligned_flux.py) | 同参数通量批次，缺模只报告部分和。 |
| [report_figure1_background.py](../src/report_figure1_background.py) | 88模场批次、断点续算、源码指纹校验。 |
| [render_figure1_background.py](../src/render_figure1_background.py) | 完整模态检查、复场合成、固定色标和原图像素对照。 |
| [snapshot_figure1_run.py](../scripts/snapshot_figure1_run.py) | 不干扰运行进程，导出可校验的阶段性发布快照。 |

详细符号、频率、角函数归一化和通量符号见 [统一约定](li_alignment/CONVENTIONS_ZH.md)。本轮为保护正在运行的源码指纹，将增补注释放在本说明中；运行模块原有 docstring 和关键公式注释保留。

## 缓存修复

`lorenz_kappa.scalar_resolvent` 容量由96增到512。L<=20 的完整扫描需要100组辅助径向解；原容量会在循环扫描时反复淘汰。修复只改缓存容量，不改数学算子。模拟扫描为100次构建、7900次复用；实际解与不经缓存的新积分逐元素一致。[修复日志](li_alignment/CACHE_CAPACITY_FIX_20260925.md)。

旧87个度规半径在验证张量数组一致后兼容迁移，[迁移清单](li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512/cache_policy_migration.json)明确标识其来源，没有冒称全部重算。

## 运行及验证

在支持 pybhpt 的 Linux/WSL Python 环境：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-environment.txt
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python src/report_li_aligned_flux.py \
  --orbit 20 --cloud-ell 1 --workers 4 --output outputs/my_li_rp20
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python src/report_figure1_background.py \
  --workers 4 --output outputs/my_figure1
```

不要在已有后台实例运行时重复启动同一输出目录。源码或参数不一致时，续算会拒绝混用数据；不要删除这个检查来强行续算。

FIG.1 自动原图比较需要本地 `outputs/` 下的原论文图像及此前校准输入，这些论文原始素材和大型缓存不随 Git 上传。其来源为 arXiv:2501.09806v1 的 `Figures/BosonEMRIsFieldPlots.pdf`，校准路径及散列见 [原图清单](environment_reproduction/figure1_pixel_comparison_20260917.json)。现有脚本按固定素材路径及散列验证；仅克隆仓库尚不能完成原图自动后处理，须先恢复这些素材。字段计算本身无需 Zotero。最终向 Windows 工作区复制结果的目标路径是当前主机专用设置，异机运行需在新批次开始前调整。

`li_source_terms.json` 已提交，通常无需重生成。若要重新生成符号系数，可将 SymPy 1.13.3 安装到独立的 `outputs/li_symbolic_dependencies/` 后运行 `scripts/generate_li_source_terms.py`，避免改变正在运行的环境。

本轮相关测试63项通过，覆盖缓存、Leaver、半解析源、BL归一化、Green核、负频率、模态完整性及复场相干叠加；14条警告为 Matplotlib/Pyparsing 弃用提示。测试通过不等于已证明全部物理/数值收敛。

## 尚未解决

- FIG.1 的88模计算与最终逐点比较尚未完成。
- 静态 completion 的有限阶跳跃误差须保留并检查，不能从通量一致推断静态场正确。
- 作者复频冻结顺序、BL质量数学视界处方、作图坐标嵌入和部分历史代码细节未完全公开。
- 度规保护阶、角向/径向积分、源外截断及 Green 边界仍需联合收敛检验。
- 尚未完成四极云全参数扫描，以及轨道反作用和累计相位演化。
