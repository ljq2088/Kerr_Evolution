# 论文与实现的18阶截断基底尚未核对

原始v1 TeX的Numerical Procedure and Validation / Metric data段说，
作者的包输出自旋加权球谐数据，然后对ell求和至18。当前
LorenzMetricMode则将nonstatic_metric的分离模态编号从|m_g|累加
至18；其中真空角函数为SpinWeightedSpheroidalHarmonic。

这说明两边具有相同的整数上限，但还不能证明有限截断算子相同。
对于分离角函数，球谐与椭球谐的关系为

\[
{}_sS_{Lm}(a\omega,\theta)
=\sum_{j\geq\max(|s|,|m|)}b^{(s)}_{jL}(a\omega){}_sY_{jm}(\theta).
\]

若某个重构分量具有展开系数A_L，则先截断分离模态给出

\[
\sum_{L\leq18}A_L\sum_j b^{(s)}_{jL}{}_sY_{jm},
\]

而先转换再截断球谐阶给出

\[
\sum_{j\leq18}{}_sY_{jm}\sum_L b^{(s)}_{jL}A_L.
\]

旋转时b并非对角矩阵，两种有限和一般不同。完整度规重构还
包含微分算子和依赖theta的几何系数，因此实际转换不能简化成
给所有BL坐标分量直接做同一个标量球谐投影。

尚待核对作者所指球谐数据的具体分量、标架、重投影位置以及
分离模态内部截断。原文这一段本身不足以唯一确定这些细节。
这项差异目前也没有数值证据能够解释无穷远5.04%、视界32.66%
或最高阶弱通道26..33%的偏差，不作此因果结论。

因此所有现有Lg=18数据继续保留，准确含义是“实现的重构模态
编号上限18”，不是“已证明与论文完全相同的度规截断”。下一步
需要分辨率和基底转换的独立比较，而不是调整通量归一化。

## 单个角函数的数值检查

report_metric_basis_mixing.py计算了rp=20M、mg=2,5,7,9,11、
s=-2..2、L=12,18的球谐展开。所有系数范数核对为1，结果见
metric_basis_mixing.json。对s=-2、L=18，j>18尾部的L2范数在
mg=5时为0.0024543733，在mg=7时为0.0033132213。这直接表明
保留L=18的椭球谐函数仍包含非零的j>18球谐成分。

该数值仅是单个归一化角函数的尾部范数，尚未经过度规重构、
径向振幅加权或源项收缩，不能解释为度规误差或通量误差。
也没有计算分离阶L>18向低球谐阶的反向贡献。

来源：outputs/environment_reference/main_PRL.tex，第681行；
实现：src/environment_lorenz_mode.py及src/lorenz_metric.py。
