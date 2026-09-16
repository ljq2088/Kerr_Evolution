# Figure 2 新增半径：r_p=15M（2026-09-16）

本次独立增加 \(\alpha=0.3\)、精确阈值 Kerr 的 \(r_p=15M\) 数据；使用当前源码重新计算完整 \(L_g=18\) 的 \(m_g=+1\) 度规，并由实度规的负 Fourier 分量关系获得 \(m_g=-1\)。没有把旧三个半径的通量插值成新结果，也没有重标定振幅。

**这里只得到 6 个实际标量通道，不是新的全模态总通量点。** 按论文使用的有限范围计数，视界 \(\ell\le5\) 覆盖 5/18，无穷远 \(\ell\le6\) 覆盖 3/9。JSON 中两端 `total_signed_flux` 均为 `null`，遗漏模式保留为 missing；作图应明确使用 partial 标记。

全部值的单位为 \(q^2(M_c/M)\)：

| 边界 | 已算部分和 | 已算/所需模式 |
|---|---:|---:|
| 视界 | -9.657087728188e-05 | 5/18 |
| 无穷远 | 5.083746471675e-07 | 3/9 |

必须区分参考版本。以下是各版本图读的**总通量幅度**，不是本地已算模式的同范围总和：

| 参考 | 视界总量幅度 | 无穷远总量 |
|---|---:|---:|
| Dyson arXiv 2501.09806v1（历史曲线） | 7.146330557308e-05 | 1.286107068550e-05 |
| Li v2 图中 later Dyson 曲线 | 9.681557947626e-05 | 1.285807880808e-05 |
| Li arXiv 2507.02045v2 曲线 | 9.249498275296e-05 | 1.332888069642e-05 |

视界部分和 / later Dyson 总量为 0.99747249，亦即部分和幅度低约 0.25275%。**这不是总通量误差**：仍有 13 个模式未算，而且此次半径的数值收敛尚未独立完成。相对 Li 总量，部分和幅度约为 1.044066 倍；这也只是 partial/total 比较。旧 v1 的比值 1.351335 仅用于历史版本对照，不能据此声称相对更新后的论文仍有 35% 偏差。

对于目标 threshold Kerr，\(F^H_{\ell m}=4Mr_+m_g^2\Omega_p(\Omega_p-\Omega_H)|Z_H|^2\le0\)，因此补算其他模式会使本地视界部分和的幅度增加；这不能替代遗漏模式的实际计算，也不能消除图读和不同参考参数/约定带来的比较限制。

无穷远部分和仅约为 later Dyson 总量的 3.95%、Li 总量的 3.81%。未算的其他 m 扇区很重要，不得把这些比值写成无穷远总通量误差。

实际通道如下；\((6,2)\) 的视界值虽也由解给出，但超出本次视界 \(\ell\le5\) 覆盖定义，因此未加入上表的 H 部分和。

| scalar (ell,m) | horizon orbital flux | infinity orbital flux |
|---|---:|---:|
| (2,2) | -1.164389887416e-06 | 4.719748762364e-07 |
| (4,2) | -2.986113020431e-09 | 3.639481305768e-08 |
| (6,2) | -1.301086027755e-11 | 4.957873460604e-12 |
| (0,0) | -9.528459474114e-05 | -0.000000000000e+00 |
| (2,0) | -1.182178757461e-07 | -0.000000000000e+00 |
| (4,0) | -6.886645628832e-10 | -0.000000000000e+00 |

径向采用 96 个源积分点：近视界对数分片 32 点，其他分片各 8 点；角向 18 点；源截断为 \(r_++0.0005M\) 到 \(320M\)，Green 外边界 \(1000M\)，视界偏移 \(10^{-4}M\)。此次新增点没有完成独立的径向/角向/L_g 全部收敛测试，因此仍标为 finite resolution、not converged。没有为这个点额外运行昂贵的绘图或补齐所有 m 扇区。

**缓存与源码谱系。** 当前源哈希为 `4b3aa366551ace4902f9d9c3669859363645ec57a14cd8650b1ddbecdd39b863`。新半径原无结果；本轮主采样器和独立后半段采样器分别生成原子缓存，在同一 metadata/source hash 验证后合并。只有完整 96 点都已核验存在，才结束仍在重复采样的旧进程，并重新加载缓存完成六个标量投影。因此最后 batch 显示 `cached_radii=96`，指本轮刚生成的数据，不是把 traceguard 变更前的历史缓存改名复用。实际过程保存在 `figure2_rp15_auxiliary_cache_20260916.json` 和 `figure2_rp15_cache_handoff_20260916.json`。

可直接读取的总表：`docs/environment_reproduction/figure2_rp15_partial_20260916.json`，含准确计数、部分和、missing modes、六通道文件列表、所有输入 SHA-256 和逐文件当前 source provenance 验证。

配套覆盖清单：`flux_coverage_L18_nt18_rp15_partial_20260916.json`。历史参考 `paper_figure2_rp15_20260916.json` 使用 `report_figure2_radii.py` 从 Dyson arXiv v1 Fig.2 矢量 PDF 独立提取。更新参考读取 `docs/paper_reproduction_20260916/figure2_published_reference_curves.json` 中的 later Dyson / Li 曲线，并在半径与 log(flux) 上插值，乘 alpha^6 转为上述归一化；来源追溯至 `li_reference_comparison_20260916.json` 所列 Li v2 原始矢量 PDF。两者均无作者原始数值表或严格图读误差界。新报告保存参考文件 SHA-256，避免混用版本。WSL 无 pdfplumber，使用已有 Windows bundled Python 执行该提取脚本；未安装新环境。

六个响应源文件均位于 `docs/environment_reproduction/`：

- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg1_sl2_inner0.0005_outer320_log_h32.json`
- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg1_sl4_inner0.0005_outer320_log_h32.json`
- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg1_sl6_inner0.0005_outer320_log_h32.json`
- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg-1_sl0_inner0.0005_outer320_log_h32.json`
- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg-1_sl2_inner0.0005_outer320_log_h32.json`
- `forced_mode_nr8_nt18_L18_alpha0.3_rp15_mg-1_sl4_inner0.0005_outer320_log_h32.json`

新增计算与报告均未修改共享生产模块，未提交 Git。
