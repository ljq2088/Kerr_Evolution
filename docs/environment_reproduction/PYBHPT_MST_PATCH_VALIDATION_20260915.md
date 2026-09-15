# pybhpt 1.0.0 低频 MST ν 分派修复与隔离验证

2026-09-15。本次确实修正了一个算法实现错误，但仅在隔离诊断库中验证；没有覆盖 `.venv`、修改其他项目或全局安装软件。它不等同于修复当前生产 AUTO 的32.66%视界通量偏差。

## 定位与补丁

检查版本是 v1.0.0，commit `9fe9c57e2d1c92d944ba70ba3c1b81b665e9d274`。原始 `nusolver.cpp` SHA256 为 `0688632c4d54e4f4faa1fc5267c7a35bce51ae524b86b6c02bfa35d5d53dbb28`。

原始代码的低频捷径存在两类错误：

1. 标量 s=0 的 ell 分派错位：ell=2、3、4 分别调用 ell=0、1、2 的专用公式；ell=0、1 则掉入在这些阶数有奇异分母的通用公式。ell=3、4 的正确专用函数已经存在，却没有调用。
2. s=±1 根本没有已实现的低频 ν 展开，但同一捷径仍调用它，并把返回的0作为 ν 使用。

在实际 a=0.8771530275949366、ω=0.011071760582072043、m=ell=1 下，原版 s0 返回 NaN，s±1 返回0。补丁修正所有 scalar ell0…6 的分派，并让 s±1 进入已有的 monodromy 初值＋MST 连分式校正路径；缺少展开的自旋不再在 monodromy 失败后的备用分支调用不存在的展开。频率范围判断使用 |ε|。此外，原 `nu_solver_guess` 对负 ε 会把正的误差估计乘成负数，导致提前接受初值；补丁也将这个判断和缩放改用 |ε|。

补丁：[`patches/pybhpt_v1_nu_low_frequency.patch`](patches/pybhpt_v1_nu_low_frequency.patch)。上游来源：[nusolver.cpp 固定提交](https://github.com/znasipak/pybhpt/blob/9fe9c57e2d1c92d944ba70ba3c1b81b665e9d274/cpp/src/nusolver.cpp)。这不是把有限 ν 无条件当作正确解：下面另行检验它满足收敛的连分式以及所产生的物理边界解。

## 实际编译与隔离范围

构建产物在 `outputs/pybhpt_mst_validation/libpatched_mst.so`。仅编译补丁后的 ν 翻译单元和一个薄 C++ 桥接层；MST 特殊函数/级数核仍链接到已安装且未经修改的扩展。桥接库使用 `-Bsymbolic` 绑定内部修正函数，用 `RTLD_LOCAL` 加载；通过非零 ν 的 `MstParameters` 显式构造原 MST workspace。

这是隔离的算法验证组件，不是完整重新编译/安装的 pybhpt Python 包。没有使用 LD_PRELOAD，没有改变生产 Python 进程的函数绑定。构建和运行前后均核对安装扩展 SHA256 未变。完整命令、构建脚本/源/头文件/补丁/二进制散列见 `outputs/pybhpt_mst_validation/build_manifest.json`，其内容同时嵌入结果 JSON。

```bash
.venv/bin/python src/build_pybhpt_mst_validation.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_pybhpt_mst_patch_validation.py
```

仅下载固定版本头文件并使用已有 g++、GSL 开发文件；未创建或安装新的 Python/Julia 环境。

## 独立 ν 与中间域物理验证

21组参数包括：Kerr s0 ell0…6，s±1 ell1、2、4，s±2 ell2，以及 Schwarzschild s0、±1 ell1 三个控制组，另外包含 Kerr s0、±1 ell1、m=-1、ω=-Ωp 三个实际负频率对照。统一采用上面的实际低频；m=min(ell,1)。s0 ell0,m0,ω≠0 是合法齐次径向方程的算法控制例，不声称它是该圆轨道实际源模态。

连分式由 Python/mpmath 重新实现，使用65位运算、双侧截断深度32/64/96，并独立求根。其公式是

\[
\beta_0+\alpha_0 R_1+\gamma_0 L_{-1}=0,\quad
R_n=-\frac{\gamma_n}{\beta_n+\alpha_n R_{n+1}},\quad
L_n=-\frac{\alpha_n}{\beta_n+\gamma_n L_{n-1}}.
\]

全部21组 ν 检查通过：补丁值与独立根的最大绝对差为 **7.50×10^-15**；补丁 ν 的独立归一化连分式残差最大 **4.20×10^-10**，这是接近抵消的连分式各项相对于其和的条件数影响。独立根在深度96处绝对残差不超过 **5.29×10^-65**，32→96的截断比较稳定。ν仍使用双精度角向本征值作为输入，此检查没有声称该输入也有65位精度。

实际 Kerr ell1 的修正结果：

| s | 修正 ν |
|---|---:|
| 0 | 0.9996911721643061 |
| -1 | 0.9996175886816967 |
| +1 | 0.9996175886816967 |

18组 Kerr（包含3组负频率对照）在 r=2.1、3、10、20、40 的 unit-transmission MST 与原 AUTO 比较：

| 指标 | 所有 Kerr 对照中的最大值 |
|---|---:|
| R 相对差 | 6.30×10^-13 |
| dR/dr 相对差 | 5.56×10^-13 |
| 加权 Wronskian 沿半径相对漂移 | 3.38×10^-13 |
| 固定 r=20 的 Green 核相对差 | 5.49×10^-13 |
| 视界辐射核 Rup(r)/W 相对差 | 5.50×10^-13 |
| 无穷远辐射核 Rin(r)/W 相对差 | 1.36×10^-13 |

这里 W=Δ^(s+1)(Rin Rup'−Rup Rin')，取稳定参考点 r=20 归一；Green 核与 AUTO 对照采用同一源测度因子。比较保留复相位，没有任意拟合整体振幅。这同时检验了原生 MST 的正确 In/Up 边界分支在中间域恢复正常。

## 仍存在的端点求值缺陷

修正 ν 后，直接 MST 的某些特殊函数/级数求值仍在端点附近失效。因此完整径向回归记录为 **18/21通过**，不是21/21：三个 Schwarzschild 控制组虽 ν 正确，但 r=2.1 接近其 r+=2 视界，Up导数仍出现 NaN 或约1.63%的误差。

实际 Kerr 的 ell=m=1，额外端点压力测试给出：

| s | r=1.5 的 Up 行为 | r=316 的 In 行为 |
|---|---|---|
| 0 | R偏约1.90%，导数NaN | R'偏约1.04% |
| -1 | R偏约7.40%，R'偏约64.7% | R/R'约1.3×10^-11 |
| +1 | R偏约1.28×10^-4，导数NaN | R偏1.27×10^-4，R'偏7.20×10^-4 |

这些是修正 ν 后仍留在原 MST 核内的另一个误差源，不能通过“强制低频全部使用 MST”解决。当前证据支持在经验证的重叠区使用 MST，并通过稳定的变换方程/ODE 继续到端点；必须重新核对 Wronskian 和源积分后才可作为生产路径。此处没有实施或声称该全域替换。

完整逐模态复值、有限性、R/R'、连分式和 Green 数据及验证脚本SHA256见 [`pybhpt_mst_patch_validation.json`](pybhpt_mst_patch_validation.json)。原版警告 `No low-frequency expansions` 出现在日志，是诊断显式调用原版 ν API 留下的对照信息；不代表隔离补丁仍走该错误分支。
