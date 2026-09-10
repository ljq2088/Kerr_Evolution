# Kerr 双曲切片标量场演化

独立 WSL 项目：`/home/ljq/code/kerr-hyperboloidal`。分支：`Kerr双曲切片演化`。

在固定 Kerr 背景上求解无质量、无源、最小耦合标量场 `Box_g Phi = 0`。
采用 horizon-penetrating hyperboloidal minimal gauge，将未来零无穷置于 `sigma=0`，未来事件视界置于 `sigma=1`。
实际演化正则场 `u=r Phi` 和 `p=partial_tau u`，包含两个端点，无人工反射边界条件。

框架入口为用户 Zotero 收藏夹“双曲框架”中的 Macedo–Zenginoglu (2025)
《Hyperboloidal Approach to Quasinormal Modes》。Kerr 系数直接核对该收藏夹中
Ripley (2022) 的式 (11)，取场自旋权重 `s=0`；这里不是将 Schwarzschild 势简单替换成 Kerr 势。
角向采用固定 `m` 的球谐 Galerkin 展开，完整保留截断空间内的 `l ↔ l±2` 耦合及所有 `im a` 项。

## 运行

```bash
cd /home/ljq/code/kerr-hyperboloidal
bash run.sh --spin 0.7 --m 2 --initial-l 2 --lmax 6 --n 64 --tmax 80 --output outputs/my_kerr
```

`--spin` 是 `a/M`，`--mass` 是 `M`，`--tmax` 是有量纲几何时间 `tau`（默认 `M=1`）。
`--n 64` 表示 65 个 Chebyshev–Lobatto 径向点；`--lmax` 必须通过收敛检查选取。
初始场是 `sigma` 上两端为零的 Gaussian 型 `l=initial_l` 脉冲，初始 `p=0`；
因此它通常含向内及向外传播两部分，不能称为纯出射脉冲。

输出包括 `evolution.npz`（时间、网格、所有球谐系数及其时间导数）、
`waveforms.csv`（两个端点各模态实部与虚部）、`waveforms.png` 和 `run.json`。
视界处的物理场是 `u/r_plus`；零无穷处记录的是辐射场 `lim r Phi`。
一个复 `m` 扇区的实部可以构造实标量场；任意三维数据可线性叠加多个 `m` 扇区。

## 验证与重现示例

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 .venv/bin/python src/validate.py
```

验证覆盖 Ripley 方程系数、球谐投影、端点因果方向、Schwarzschild 解耦、
Kerr 模态耦合、`m ↔ -m` 共轭关系及质量缩放。
收敛脚本比较径向 `N=32,48,64`、角向 `lmax=4,6,8` 和时间容差，并生成 `m=0,2` 两组 Kerr 示例。
定量结果见 [验证记录](docs/validation.md)。

## 环境

项目使用独立 `.venv`，不共享或修改其他项目的第三方包。重建：

```bash
/home/ljq/miniconda3/envs/few_env/bin/python -m venv .venv
.venv/bin/python -m pip install -r requirements-lock.txt
```

源文献不需要随运行加载。依赖已锁定，运行不需要 Zotero 或网络。
详细推导见 [Kerr 方程说明](docs/kerr_hyperboloidal_zh.md)，文献阅读范围见
[Zotero 来源记录](docs/zotero_sources.md)。

## 当前范围

实现针对 `|a/M|<1` 的线性测试场。没有加入标量质量、源项、自相互作用、背景反作用或 EMRI 粒子源。
本次验证针对给定初值、`a/M=0.7`、`tau<=80M` 的波形收敛；并不自动证明任意参数或极长时间尾波的精度。
图中的衰减振荡是时域结果，本项目尚未实施 QNM 频率拟合、频谱识别或能量/角动量通量守恒检查。
