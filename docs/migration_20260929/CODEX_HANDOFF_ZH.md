# Kerr / EMRI 环境效应项目：新电脑 Codex 交接手册

交接日期：2026-09-29。本文可以单独发送给新电脑上的 Codex；无需旧对话上下文。

## 1. 给接手 Codex 的任务

请先阅读本手册，恢复项目、数据和 Python 环境，完成最低限度验收，然后汇报实际恢复情况和下一步建议。遇到本机已有项目时先检查，不能覆盖已有修改。

项目长期目标是严格按照论文对象、参数、方法、符号与归一化，复现 Kerr 黑洞周围标量云受圆轨道伴星扰动的场和通量；以 Li 后续修正为主要基准，同时区分 Dyson 原图。当前阶段是迁移与接手，不自动重启已取消的长时间扫描。

用户工作习惯：

- 用中文交流；原理解释需要详细、严谨，明确代码对应与证据边界。
- 已授权的常规操作直接完成，不反复询问确认。
- 长计算放后台，写检查点和日志；不要在对话里不断轮询。若只剩等待，结束当轮并让后台继续，不持续消耗 token。
- 不擅自增加定时提醒、自动化或子代理。
- 不通过任意乘倍数、调相位、换归一化或省略难算模式来宣称复现成功。
- 旧任务存在取消/停止记录；不要因看到历史 `running` 或脚本里的自动串行计划就启动它。

## 2. 唯一明确的恢复来源

| 内容 | 地址或标识 |
|---|---|
| GitHub 仓库 | https://github.com/ljq2088/Kerr_Evolution |
| 当前工作分支 | `codex/reproduce-kerr-environment` |
| 迁移快照标签 | `migration-2026-09-29` |
| 标签对应的代码提交 | `639213733d9a0da1351431df7e17cab4bdb4a75d` |
| 上传回执提交 | `6c705df0aea312a1341096421544241f68b20a98`，后续可能还有文档提交 |
| 完整数据下载页 | https://github.com/ljq2088/Kerr_Evolution/releases/tag/migration-2026-09-29 |
| 仓库内交接记录 | `docs/migration_20260929/` |

**默认分支 `Kerr双曲切片演化` 不是当前环境效应工作分支。只克隆默认分支或只下载源码 ZIP 不够。**

迁移包共约783 MB：

| 文件 | 字节数 | 内容 |
|---|---:|---|
| `kerr-project-state-20260929.tar.gz` | 539691510 | 31,468 个项目文件：完整工作树、原先忽略的 outputs、数值结果、检查点、参考素材、度规缓存、后续实验脚本 |
| `kerr-workspace-reports-20260929.tar.gz` | 242893787 | 876 个 Windows 工作区文件：中文报告、推导、图、审查和临时实验 |
| `SHA256SUMS` | 206 | 上述两个包的整体 SHA256 |
| `assets.json` | 1949 | 文件大小、散列、数量与排除记录 |

项目包 SHA256：

```text
5f71741b0c0f08055f0abc48a41e1020b94982981fd12f0eafdf4dfca13712f0
```

工作区报告包 SHA256：

```text
3247e43747215977b489ac3ed7600bc64935dec9952e356d61d96cdb150c7ce0
```

GitHub 返回的大小与 SHA256 已在原机核对。每个包内还有逐文件 `MIGRATION_FILES.json`，恢复工具会验证它。

未打包 `.git`、`.venv`、字节码、测试缓存、egg-info 和两个可重装依赖目录；Git 历史由 clone 获取，Python 环境重新安装。系统程序、Wolfram 许可证、账户凭据、Zotero 整个数据库和运行中的内存不在备份中。原作者公开参考代码及已有论文素材保留各自许可。

## 3. 新电脑恢复步骤

### 3.1 检查环境

推荐 Windows + WSL2 Ubuntu x86_64，或原生 Linux x86_64。原机 Python 为 **3.12.2**，核心依赖为 numpy 1.26.4、scipy 1.14.1、pybhpt 1.0.0、mpmath 1.3.0。完整版本见 `requirements-migration.txt`。

先检查 Git、curl、Python 3.12 和 venv 支持。若系统的 `python3` 不是3.12，请安装/选择3.12，不要直接升级系统 Python，也不要复制旧虚拟环境。若新机为 ARM/macOS，先核实 pybhpt 的安装与数值兼容性，不默认它等价于原 x86_64 环境。

建议放在 Linux 主目录 `~/code`，避免把主计算仓库放在 WSL 的 `/mnt/c` 文件系统下。

### 3.2 克隆与下载

下面命令仅适用于目标目录尚不存在的情况；已有目录先检查其分支、修改和文件，不可直接覆盖。

```bash
mkdir -p ~/code
cd ~/code
git clone --branch codex/reproduce-kerr-environment https://github.com/ljq2088/Kerr_Evolution.git kerr-hyperboloidal
cd kerr-hyperboloidal
git switch --detach migration-2026-09-29
git rev-parse HEAD
```

确认 HEAD 为 `639213733d9a0da1351431df7e17cab4bdb4a75d`。暂时固定迁移标签是为了避免将旧数据恢复到未来已经变化的代码上。

```bash
mkdir -p ../kerr-migration-download
cd ../kerr-migration-download
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/kerr-project-state-20260929.tar.gz
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/kerr-workspace-reports-20260929.tar.gz
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/SHA256SUMS
curl -fLO https://github.com/ljq2088/Kerr_Evolution/releases/download/migration-2026-09-29/assets.json
sha256sum -c SHA256SUMS
```

必须两项都为 OK；下载失败时不要把 HTML 错误页当作压缩包。

### 3.3 验证并恢复

```bash
cd ../kerr-hyperboloidal
python3.12 scripts/restore_migration_archive.py ../kerr-migration-download/kerr-project-state-20260929.tar.gz --destination .
python3.12 scripts/restore_migration_archive.py ../kerr-migration-download/kerr-workspace-reports-20260929.tar.gz --destination ../EMRI-workspace-reports
```

工具先校验所有文件，然后只写入缺失文件；已有内容相同则跳过，已有内容不同则拒绝。遇到冲突先调查，或换一个干净克隆，不要绕过保护。若只验证而暂不恢复，省略 `--destination`。

两个包解压总内容约1 GB，还需为 Python 环境、Git 和后续缓存预留空间。工作区报告与项目分开恢复，避免混淆同名 `outputs`。

### 3.4 重建依赖

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-migration.txt
```

锁定清单包含本项目两个本地可编辑包，必须在项目根目录安装。不要未经记录改成最新依赖。安装失败时记录 Python 版本、平台和具体依赖错误，再针对问题处理。

仅在需要重跑符号生成/HDF5原作者审查时安装：

```bash
.venv/bin/python -m pip install --target outputs/li_symbolic_dependencies sympy==1.13.3 mpmath==1.3.0
.venv/bin/python -m pip install --target outputs/paper_metric_reference/python_deps h5py==3.16.0 --no-deps
```

主 Python 场/通量计算不要求 Zotero 或 Wolfram 实时运行。重跑 Mathematica 原作者程序则需要新电脑自己的安装和授权；已导出的数据在迁移包中。

### 3.5 最小验收

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_figure1_background.py tests/test_li_flux_pipeline.py tests/test_li_alignment_contract.py tests/test_li_flux_inventory.py
```

原机结果：**15 passed**，14条 Matplotlib/Pyparsing 弃用提示。测试通过不代表论文已复现。

核对完整88模和源码指纹；以下脚本只读取文件，不重新计算：

```bash
.venv/bin/python - <<'PY'
import hashlib,json
from pathlib import Path
root=Path.cwd()
d=root/'docs/figure1_alignment_20260926'
e=json.loads((d/'execution_completed_20260929.json').read_text())
expected={(l,m) for l in range(2,13) for m in range(-l,l+1) if (l+m)%2==0}
seen=set()
for f in d.glob('mode_l*_m*.json'):
    v=json.loads(f.read_text()); key=(v['ell'],v['m'])
    assert key not in seen
    seen.add(key)
    assert v['config']==e['config']
    assert v['implementation_sha256']==e['implementation_sha256']
    assert hashlib.sha256((d/v['radial_file']).read_bytes()).hexdigest()==v['radial_sha256']
    bank=d/f"metric_mg{v['m']-1}.npz"
    assert hashlib.sha256(bank.read_bytes()).hexdigest()==v['metric_bank_sha256']
assert seen==expected and len(seen)==88
for name,digest in e['implementation_sha256'].items():
    f=root/'src'/(name if name.endswith('.json') else name+'.py')
    assert hashlib.sha256(f.read_bytes()).hexdigest()==digest,name
assert (d/'figure1_fields.npz').is_file()
print('88 modes, metric/radial checksums and original implementation verified')
PY
```

验收后查看分支相对于迁移标签的变化：

```bash
git diff --stat migration-2026-09-29 origin/codex/reproduce-kerr-environment
```

若只有交接文档/回执等新增，可以切回 `codex/reproduce-kerr-environment`。若有科学源码变化，先核对新实现和缓存兼容性，必要时基于迁移标签新建本地工作分支。不要长期在 detached HEAD 上直接开发。

## 4. 项目究竟在计算什么

项目包含两条不同工作线：

1. **早期双曲切片时域演化**：Kerr 上无质量标量场，以及后续导入的 s=-2 点粒子 Teukolsky 演化。README 中部分旧图、`run.sh`、`run_sminus2.sh` 属于这一条。
2. **当前主要工作：环境效应频域计算**。固定 Kerr 背景和圆形赤道轨道，先求标量云准束缚态，再求伴星引起的 Lorenz 度规扰动，构造受迫标量源，解径向方程并相干合成场，最后计算有效轨道能量通量。

当前论文场图不是把旧的无源高斯脉冲时域演化改名。当前场图采用固定轨道，不把辐射通量反馈为 r0(t)。轨道能量平衡和累计相位是另一阶段。

参考论文：

- Dyson 等，*Environmental effects in extreme mass ratio inspirals: perturbations to the environment in Kerr*，arXiv:2501.09806v1。项目源文件：`outputs/environment_reference/main_PRL.tex`。
- Li 等，*Extreme mass-ratio inspiral within an ultralight scalar cloud I. Scalar radiation*，arXiv:2507.02045v2。源文件：`outputs/paper_original_reference/li_2507_02045v2/source/main.tex`。
- Lorenz 重构参考包括 *Metric perturbations of Kerr spacetime in Lorenz gauge: circular equatorial orbits*（DOI 10.1088/1361-6382/ad52e3）与 *Sourced metric perturbations of Kerr spacetime in Lorenz gauge*。已有公开原作者代码、笔记和对照位于 `outputs/paper_metric_reference/` 及相关 docs。

2026-09-15 已找到并运行公开作者 Lorenz 代码；不能继续沿用更早“完全没有作者代码”的说法。但尚没有该环境场图在全部目标参数下的作者原始数组与完整出图程序。

## 5. 已完成结果与问题清单

### 5.1 Dyson FIG.1 场景：完成计算，尚未复现

目录：`docs/figure1_alignment_20260926/`。

- M=1，a=.88，mu=.3，r0=3.5；背景云 (ell_c,m_c,n_c)=(1,1,0)，即 |211>。
- 标量 ell=2..12，共88个允许模态；包含正负方位模态、束缚通道和静态度规驱动通道。这里“静态驱动”不等于标量总频率为零。
- 度规椭球输入L20，加权球谐输出j18；度规角积分q40，源q64，Fourier pmax12；径向常规段24、近视界段64。
- 相干求和复场后取绝对值，未拟合物理结果的幅度或相位。
- 图片：`figure1_comparison.png`；数据：`figure1_comparison.json`；95MB复场数组：`figure1_fields.npz`；最终记录：`execution_completed_20260929.json`。

**目标和方法必须分开表述：**本次采用 Li 修正后的实现、精确a=.88，计算 Dyson 的 r0=3.5 场景，再与 Dyson 旧图对照。Dyson 图注为a=.88，而同步云背景约a=.877153；作者复频处方、静态completion和图坐标尚未全部核实。因此不是严格同处方复现，也不是 Li r0=41.1/42.1M 的场图。

对8个赤道/8个子午面色标采样点做单一倍数的最小二乘诊断：

| 采样范围 | 最佳乘数 | 乘后相对L2误差 |
|---|---:|---:|
| 赤道面 | 5.1577856 | 4.836% |
| 子午面 | 5.2359234 | 42.038% |
| 两者合计 | 5.1764094 | 22.666% |

该诊断使用图像读数，不是作者数组，也不是全图误差。一个公共倍数不能解释子午面形状差异。ell8→12的复场网格变化约1%，但不能据此宣称度规、径向、源外截断等已经收敛。不能擅自把全部结果乘5.2。

### 5.2 Li 同参数总通量

r0=20M、a=.88、mu=.3、|211>、18个非静态通道；目录：`docs/li_alignment/rp20_cloud11_j18_L20_q40_nr12_h64_cache512/`。

以下单位为 q²(Mc/M)，q为粒子/黑洞质量比：

| 边界 | 本地有效能量通量 | 对Li图读数的差异 |
|---|---:|---:|
| 无穷远 | 1.8978444432481997e-5 | +0.186696% |
| 视界 | -6.828254573842368e-5 | +3.220241% |

另一个请求 r0=18.3M 已于2026-09-29完成18模，记录在 `outputs/li_fig2_20260928/r18p3_completion.json`。**该文件单位是 q² zeta²，与上表相差alpha^6，不能直接混用：**无穷远0.0235529372566338，视界−0.11228157696378695；与对应图读数相差约+0.1663%/+3.7469%。用户只要求完成此半径后停止，未完成全半径扫描。

Li 指出 Dyson 的标量分离常数和视界通量归一化存在修正。早期相对旧视界图的约32.7%差异，不能再不加区分地称为当前实现未解释的误差。

### 5.3 其他后续工作

| 目录 | 已记录状态/限制 |
|---|---|
| `outputs/li_fig45_20260927/` | FIG.4/5 源项有限分辨率比较已完成，查看 comparison 与脚本具体定义 |
| `outputs/li_fig9_10_20260928/` | 部分逐模结果与归一化诊断；剩余扫描用户已取消 |
| `outputs/dyson_fig7_20260927/` | 粒子处多极图任务 incomplete_with_failures；不是Li场图FIG.7 |
| `outputs/spectral_radial_trial_20260928/` | 相同固定源下，径向谱BVP与Green通量约1e-11量级吻合；未验证度规正确 |
| `outputs/li_dipole_check_20260928/` | 保留诊断数据；execution历史running但原PID已不存在，不能认为任务在运行或已完成 |

Li FIG.10 的归一化诊断：未平移、a0=1的原始 Leaver 级数相对于质量归一化云，预测通量因子 K_211≈463.9105633、K_322≈1.0541490178e8。该因子来自背景云，不是拟合图线；能解释部分无穷远通道的显示尺度，但未解释所有视界模式，也未证实作者实际出图约定。详细记录 `FIG10_FOLLOWUP.md`。不要把这个通量因子与FIG.1的场幅倍数5.2混为一谈。

## 6. 统一物理约定

先读 `docs/li_alignment/CONVENTIONS_ZH.md`，不要仅凭符号名称推测定义。

- G=c=M=1；alpha=mu M；q=mp/M；zeta=alpha³ sqrt(Mc/M)。不同论文中 epsilon 的用途不同，转换时先写清物理意义。
- 复标量场，角函数采用全立体角归一化为1的约定。
- |211> 空间云使用复Leaver频率，时间频率取其实部；本批 omega_c≈0.29629353472811226+2.216612802247371e-9 i。
- 云采用有限 Boyer–Lindquist t=0 Killing 能量归一化；质量积分近视界截断1e-6，外端2000M。作者数学视界处方尚未完全确认。
- 驱动频率 omega=Re(omega_c)+(m-m_c)Omega_g。
- 有效轨道能量通量包含云的电荷/质量记账：F=(omega-Re(omega_c))N，不能与普通场能量omega N混淆。
- 单位云质量场转换到论文系数时除alpha³；单位云质量通量转到q² zeta²单位时除alpha^6。
- 在需要求和的场图中先加复模态，再取模。纯通量比较通常不能验证相对相位、静态通道或局域场正确。
- 源项实现补回协变恒等式要求的混合迹项；Li印刷附录控制开关用于诊断。印刷式问题不证明作者实际程序使用了错误公式。

## 7. 代码导航

| 模块 | 职责 |
|---|---|
| `src/li_leaver_cloud.py` | 准束缚复谱和Leaver径向级数 |
| `src/li_normalized_cloud.py` | 云质量归一化、复空间剖面/实时间频率处方 |
| `src/environment_lorenz_mode.py` | 本地Lorenz度规重构 |
| `src/environment_static_lorenz.py` | 静态匹配与completion |
| `src/environment_metric_sampling.py` | 度规采样与缓存 |
| `src/li_separated_source.py`、`src/li_source_terms.json` | 半解析源及符号系数 |
| `src/li_order4_green.py` | 正频率四阶边界Green求解 |
| `src/li_field_green.py` | 带符号频率的延迟边界处理 |
| `src/environment_response.py` | 连续径向场重构 |
| `src/report_li_aligned_flux.py` | 同参数通量批次 |
| `src/report_figure1_background.py` | 88模场批次与完整性校验 |
| `src/render_figure1_background.py` | 复场合成、固定色标、像素比较 |
| `src/source_provenance.py` | 缓存/检查点实现指纹 |
| `scripts/restore_migration_archive.py` | 数据包校验与保守恢复 |

后续实验驱动也可能在 `outputs/` 中，不能当垃圾清理。例如 `outputs/li_fig2_20260928/run_only_r18p3.py`、`outputs/li_fig9_10_20260928/run_comparison.py`、`outputs/spectral_radial_trial_20260928/trial.py`。`run_priority_sequence.py`会串联大任务，当前不应执行。

## 8. 迁移陷阱与继续工作的原则

1. 原机路径 `/home/ljq/code/kerr-hyperboloidal` 与 `C:\Users\赖景祺\Documents\ChatGPT\EMRI环境效应` 仅供识别来源。新机使用实际路径，不假设用户名相同。
2. 历史JSON中绝对路径和PID是来源记录；不要批量替换以“修复路径”，会改变文件散列。必要时新增路径映射/新运行记录。
3. 旧FIG.1渲染器结尾硬编码Windows导出路径并使用copy2。原机曾在时间戳复制时报错，图和数组早已算完；已恢复导出。看已有图无需重新渲染，更不需要重算88模。若以后参数化导出，应建立新代码版本和新输出目录，不伪装同实现续算。
4. 缓存容量96→512是已完成的性能修复，不改变算子；不要退回旧容量。相关记录 `docs/li_alignment/CACHE_CAPACITY_FIX_20260925.md`。
5. 原有数据用实现散列防止混用。改动源码/JSON公式后不能删除检查来强行续算；旧结果仍可读取，新结果必须有明确来源。
6. 截止迁移未发现运行中的项目计算进程。部分README和execution是历史状态，优先看 `environment_and_batches.json` 和用户取消记录。
7. 先把论文的对象、参数、角向模态集合、显示量、云归一化、时间/相位切面及坐标映射固定成比较合同，再评估偏差。优先沿Li修正后的相同场图设置继续，而非拿不同处方图做“严格误差”。
8. 优先复用已验证数据做便宜诊断；昂贵重算前先有能区分假设的检验。静态completion残差、相对模态相位、源/度规离散与参考归一化仍是开放问题，原因未锁定。
9. 后续变更必须更新日志和参数/版本来源；区分“计算完成”“离散收敛”“与图一致”“物理公式验证”四种不同结论。

## 9. 接手后第一次汇报的内容

请完成恢复和最低验收后，向用户简要说明：

- 新机项目绝对路径、分支/提交、Python与核心库版本；
- 两个下载包校验和恢复是否成功；
- 15项测试及88模完整性/指纹核对结果；
- 能否打开现有FIG.1图片和读取完整场数组；
- 遇到的实际兼容问题，以及建议优先调查的一个具体问题。

此时不自动启动全半径扫描、不恢复已取消任务、不宣称论文已经复现。后续根据用户的新指令推进。
