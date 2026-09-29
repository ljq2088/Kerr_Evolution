# 新电脑交接与完整备份（2026-09-29）

## 先读这里

工作分支：`codex/reproduce-kerr-environment`。默认分支 `Kerr双曲切片演化` 不是最新环境效应工作，克隆时必须指定工作分支。

[迁移版本与数据下载](https://github.com/ljq2088/Kerr_Evolution/releases/tag/migration-2026-09-29)。这是迁移快照，不是科学验收通过的版本。

本次检查未发现仍在运行的项目计算进程；只保存文件与检查点，不迁移操作系统进程。不要根据历史 PID 或旧文档里的“running”自动续算。

## 科学进度与未解决问题

- **Dyson FIG.1 场景**：88/88 模态完成，r0=3.5M、a=.88、mu M=.3、偶极基态云、ell=2..12。使用 Li 修正后的实现，与 Dyson 旧版原图比较；不是 Li 的 r0=41.1/42.1M 场图，也未达到严格同处方复现。Dyson 图注 a=.88 与方法的同步云约 a=.877153 有区别。作者复频处理、质量积分视界处方、静态 completion、显示坐标尚未逐项确认。
- FIG.1 赤道/子午面采样点最佳诊断倍数约5.158/5.236，缩放后采样相对 L2 误差4.84%/42.04%。**没有把拟合因子乘进任何物理结果**。幅度和形状偏差仍在，不能声称复现成功。ell8→12 复场网格变化约1%，不等于全部空间收敛。
- **Li r0=20M 总通量**：18个非静态通道完成。无穷远1.8978444432481997e-5，视界−6.828254573842368e-5，单位为质量比平方乘 Mc/M。相对 Li 读图约+0.1867%/+3.2202%。有限分辨率与读图误差限制保留。
- **后续 Li FIG.2**：`outputs/li_fig2_20260928/r18p3_completion.json` 保存 r0=18.3M 的已完成请求；完整半径扫描未做完，用户要求仅完成这个半径后停止。
- **Li FIG.9/10**：`outputs/li_fig9_10_20260928` 保留部分逐模数据、图和诊断。剩余扫描已由用户取消，不自动恢复。原始 Leaver 级数归一化解释部分无穷远曲线倍数，是绘图约定假设，不是已证实的作者程序错误；视界通道仍有问题。见该目录 `FIG10_FOLLOWUP.md`。
- **Li FIG.4/5 源项比较**：`outputs/li_fig45_20260927` 保存已完成的有限分辨率比较。
- **Dyson 粒子处 FIG.7**：`outputs/dyson_fig7_20260927` 的状态是 incomplete_with_failures，不能与 Li 场分布 FIG.7 混淆。
- **谱方法试验**：`outputs/spectral_radial_trial_20260928` 对相同固定源，标量径向谱 BVP 与 Green 解吻合约1e-11量级；未修复视界通量差异，也未验证度规正确。
- **偶极诊断**：`outputs/li_dipole_check_20260928/execution.json` 仍有历史 running 标签，但该 PID 已不在运行。保留现有结果，状态未知部分不能当作完成。

上述目录有一些早期 README 仍使用当时的运行时态；当前状态以本交接文件和 `environment_and_batches.json` 的迁移观察为准。不改写历史日志。

## 数据在哪里

Git 保存源码、文档、所有88模径向解、度规银行、对照图片与 JSON。额外 Release 附件：

1. `kerr-project-state-20260929.tar.gz`：整个项目工作树（含原先忽略的 outputs、后续临时驱动脚本、参考素材、原作者公开代码与许可证、缓存、原始结果、完整95 MB场数组、执行状态与日志）。保留原相对路径。
2. `kerr-workspace-reports-20260929.tar.gz`：Windows 当前工作区的全部项目报告、推导、图片、诊断副本、临时实验脚本与已有压缩包。单独恢复，避免与 Linux 项目路径混合。
3. `assets.json` 和 `SHA256SUMS`：下载校验。每个压缩包内部另有逐文件 `MIGRATION_FILES.json`。

排除 .git、.venv、__pycache__、.pytest_cache、egg-info、pyc，以及可重装的 `outputs/li_symbolic_dependencies` / `outputs/paper_metric_reference/python_deps`。逐项排除清单在包内；Python环境以锁定清单重建，Git历史在远端。未打包系统安装、Wolfram许可证、Zotero整个数据库、账户凭据或本项目之外的目录。已有参考论文摘录和项目书目已保存。外部参考文件仍适用各自许可，不因迁移归本项目所有。

## 新电脑恢复

建议 Ubuntu/WSL2 x86_64、Python 3.12（原机3.12.2）；不要复制旧虚拟环境。以下在 Linux 终端运行，先安装 Git、curl 及 Python venv 支持。

```bash
mkdir -p ~/code
cd ~/code
git clone --branch codex/reproduce-kerr-environment https://github.com/ljq2088/Kerr_Evolution.git kerr-hyperboloidal
cd kerr-hyperboloidal
git switch --detach migration-2026-09-29
mkdir -p ../kerr-migration-download
cd ../kerr-migration-download
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/kerr-project-state-20260929.tar.gz
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/kerr-workspace-reports-20260929.tar.gz
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/SHA256SUMS
sha256sum -c SHA256SUMS
cd ../kerr-hyperboloidal
python3 scripts/restore_migration_archive.py ../kerr-migration-download/kerr-project-state-20260929.tar.gz --destination .
python3 scripts/restore_migration_archive.py ../kerr-migration-download/kerr-workspace-reports-20260929.tar.gz --destination ../EMRI-workspace-reports
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-migration.txt
```

恢复工具先验证所有文件，拒绝覆盖内容不同的现有文件。若有冲突，使用新的克隆/空目录，不要强行覆盖新工作。压缩包不含 Git 历史；克隆保留历史。恢复后回工作分支：`git switch codex/reproduce-kerr-environment`（如果远端以后继续更新，应先基于迁移标签建本地分支）。

可选外部验证依赖（仅需要重跑相关旧审查时安装）：

```bash
.venv/bin/python -m pip install --target outputs/li_symbolic_dependencies sympy==1.13.3 mpmath==1.3.0
.venv/bin/python -m pip install --target outputs/paper_metric_reference/python_deps h5py==3.16.0 --no-deps
```

Wolfram 原作者对照需要在新机器独立安装并授权 Mathematica/Wolfram；其源码/现有导出数据已保存。主 Python 通量和场求解不依赖 Zotero 桌面或 Wolfram 实时运行。

## 验证和继续工作的入口

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python -m pytest -q  tests/test_figure1_background.py tests/test_li_flux_pipeline.py  tests/test_li_alignment_contract.py tests/test_li_flux_inventory.py
```

主模块对照见 `docs/PROJECT_STATUS_20260927.md`，其9月27日的70/88进度已经被本次最终88模状态取代。符号/通量/归一化合同为 `docs/li_alignment/CONVENTIONS_ZH.md`。所有通量“百分比接近”结论仅适用于记录的参数、模态和参考图。

- 主频域入口：`src/report_li_aligned_flux.py`。
- 已完成的场批次：`docs/figure1_alignment_20260926`；原生成器 `src/report_figure1_background.py`。
- Li 后续驱动保留在各自 `outputs/li_*` 目录，迁移包恢复后可用。部分脚本没有通用CLI，先读代码/现有执行记录。
- 科学生产源码未因本次迁移修改；检查点的实现SHA256仍可审查。不同源码/参数不允许强行复用缓存，不可删除指纹检查。
- 历史 JSON 中的绝对路径和 PID 是来源记录，不是新机有效地址。查找文件时以归档内相对目录为准。不要批量改写记录，否则会破坏散列证据。
- **旧 FIG.1 渲染器结尾仍硬编码原 Windows 导出路径，并使用 copy2**。图和数组在这个导出步骤之前已经写好；本次原机曾仅在复制时间戳失败，已通过内容复制恢复。查看已完成结果无需运行它。新批次若要重新生成，应先把导出路径参数化并使用新输出目录；不要改完代码后冒充原88模批次的同实现续算。
- 下一阶段先严格区分 Dyson r3.5 场景与 Li r41.1/r42.1 场景，核对相同参数、云质量与图示量后再比较。当前迁移不会重启已取消的扫描。

本次迁移验证记录另见 `validation.txt`、`archive_roundtrip.json` 和 Release 校验清单。下载包备份的是截至本次归档的文件状态，不能恢复未写盘的内存计算。
