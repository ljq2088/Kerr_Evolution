# s=0,±1 的独立 MST 锚点径向审查

2026-09-15。审查对象是实际 alpha=.3、a=.8771530275949366、rp=20、ell=1、m=±1、omega=±Omega_p 的无质量度规/规范势径向函数。不是有质量环境标量的径向方程。

## 独立参考算法

使用前一轮隔离修复的 MST ν 选择器和原 MST 级数核，在其可靠的中间区域直接计算单位 Teukolsky transmission 归一化的 R、dR/dr，再由独立 SciPy DOP853 积分到实际源网格。没有调用 NativeRadial.solve、HBL/TEUK 的渐近边界初始化；NativeRadial 构造器仅提供角向本征值。

为避免 Delta=r²−2r+a² 在视界附近相减损失，使用

\[
q=r-r_+,\quad x=\log q,\quad \Delta=q(q+r_+-r_-),\quad P=\partial_xR=q\partial_rR.
\]

若径向方程写成 R''=−A R'+B R，则实际积分的是

\[
\partial_xR=P,\qquad \partial_xP=(1-qA)P+q^2BR,
\]

其中

\[
A=\frac{2(s+1)(r-1)}{\Delta},\qquad
B=\frac{\lambda-4is\omega r}{\Delta}
 -\frac{K^2-2is(r-1)K}{\Delta^2},\qquad K=(r^2+a^2)\omega-am.
\]

对每个边界解按锚点状态范数缩放 ODE 初值，返回时乘回同一范数；没有对 AUTO 拟合任何复振幅或相位。MST 锚点已确定 In/Up 物理边界支，因此 ODE 延拓只传播同一齐次解。

## 传播稳定性与锚点选择

不能盲目从同一个远处锚点向两个方向积分。初始测试使用共同锚点10和20、容差3e−11/3e−13/3e−14。s=−1 In 从20向视界传播的误差随收紧容差不单调，约停留在1e−8；从10传播约1e−9。该病态传播诊断被保留在 `gauge_radial_anchored_initial_audit.json`，没有删除或重标成最终可靠结果。

最终采用：

- In 在 **r=3** 锚定，独立换为 **r=2.1** 检查。
- Up 在 **r=20** 锚定，独立换为 **r=40** 检查。
- DOP853 相对容差 **3e−14**，对锚点归一状态的绝对容差 **3e−19**；另外比较3e−13、3e−11。
- 缓存 ODE 范围 r∈[r_++1e−5,400]；实际验证为原有全部88源节点，范围 **1.4807164758793971…316.8231885198029**，另加精确r=10、20作锚点/拼接核验。

s=−1 In 的两锚点差异降至约8e−12，表明下移锚点确实改善了传播条件性。s=+1 Up 的远区仍有约5.33e−11的两锚点差异，应将其计入参考误差，不宣称所有节点均有14位有效数字。

## 结果

全部 s=0,±1 × m=±1 × In/Up 均有限。以下是在所有实际源节点和正负频率上的最大相对量：

| 检验量 | 最终参考对 AUTO | 两套独立锚点之间 |
|---|---:|---:|
| 所有 R/R' | 2.691e−10 | 5.328e−11 |
| 视界辐射核 Rup/W | 2.690e−10 | 5.319e−11 |
| 无穷远辐射核 Rin/W | 6.661e−12 | 7.809e−12 |
| 固定观察半径r=20的分段 Green 核 | 1.308e−10 | 5.319e−11 |

最终参考加权 Wronskian

\[
W=\Delta^{s+1}(R^{\rm In}R^{{\rm Up}\prime}-R^{\rm Up}R^{{\rm In}\prime})
\]

沿半径的最大相对漂移 **3.016e−12**，r=20处与AUTO的绝对归一化相对差 **1.313e−13**。容差3e−13→3e−14导致的最大R/R'变化约 **8.99e−12**。

分段 Green 函数采用 Rin(r_<)Rup(r_>)/W，波函数连续，单位delta源的加权导数跳跃恒等式在2.31e−16以内；这是拼接代数一致性检查，不能替代上面的锚点、容差和独立初始化检查。

这些数据强烈限制了 s=0,±1、ell=1 径向函数本身的误差幅度，但是否足以排除32.66%的最终通量偏差，还需要将该参考注入完整源并传播到ZH。当前报告不提前宣称完成该最后一步。

## 注入接口与复现

```python
from environment_mst_anchored_radial import MSTAnchoredRadial
radial = MSTAnchoredRadial(s, ell, m, a, omega, radii,
                          anchor=3.0, up_anchor=20.0, rtol=3e-14)
radial.solve(bc=None)
Rin = radial.radialsolutions('In')
Rin_prime = radial.radialderivatives('In')
```

`radialsolution`、`radialderivative`、`radialderivative2`、数组版本及 `.eigenvalue` 均为NativeRadial兼容子集。`.diagnostics`仅含可严格JSON序列化的原生数据；`.provenance`带实现散列。三组分量配置应使用新进程避免源幅缓存交叉。

```bash
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_gauge_radial_anchored_audit.py
```

最终原始数据：[`gauge_radial_anchored_audit.json`](gauge_radial_anchored_audit.json)。实现 SHA256：`3b372683cf4171b92e33d73cda00a0278e8b6a732f5f337ba663be3433000f9b`。

本次新增参考模块和诊断报告；没有覆盖`.venv`、修改生产模块绑定、修改其他项目或提交代码。
