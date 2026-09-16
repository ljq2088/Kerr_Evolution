# 固定 Kerr 自旋 a=0.88 的复云谱：独立射击核验

对 mu=0.3、ell=m=1、最低径向态，独立双端射击给出

\[
M\omega=0.29629353472561143+2.216609396778975\times10^{-9}i.
\]

时间约定 exp(−i omega t) 下这是超辐射增长态，振幅增长时间约4.51e8 M。该自旋的 Omega_H=0.2983104071127746，明显不同于 Re omega；因此它不是严格同步的 stationary cloud。

射击使用完整复角矩阵的一般 eig：

\[
[\ell(\ell+1)\delta_{\ell\ell'}-c^2(\cos^2\theta)_{\ell\ell'}]b_{\ell'}=A b_\ell,
\qquad c^2=a^2(\omega^2-\mu^2).
\]

不能对复 c² 使用 eigvalsh。精化根对应 A=2.000342293290726−2.0342423984925983e−10 i。

径向方程为

\[
(\Delta R')'+\left[\frac{((r^2+a^2)\omega-am)^2}{\Delta}-\mu^2r^2-a^2\omega^2+2am\omega-A\right]R=0.
\]

左端用 ingoing Frobenius 解，右端用无穷远衰减的逆幂展开，以 DOP853 双端积分并令匹配点的 R'/R 相等。没有使用 Leaver 根作初猜；初猜来自已有同步云附近。

两档同时改变：外界1000→1555.56、视界偏移1e−5→1e−6、rtol2e−11→2e−13、Frobenius阶4→6、远场阶6→8、角基20→28、匹配半径22.222→25。实部变化2.22e−16、虚部变化1.70e−17，最终匹配残差约8.47e−16。

双方独立求根后交换结果，与另一代理的 Leaver N500 对照，按记录的十进制数，实部差约−3.28e−17、虚部差约+3.33e−19。两类独立算法一致；连分式自身截断变化另存于对方报告。这不是仅把求根器的内部残差当成误差估计。

仅按 Re omega 推出的 m=2 电离半径为41.6627062842 M。没有重算云质量归一化、度规、协变源或强迫标量场。若后续引入该复云，必须明确视界穿越切片上的归一化和是否冻结 Im omega；不能直接套用同步云的规则视界初值及其 BL 归一化积分。

复现脚本：`src/report_fixed_spin_kerr_cloud_shooting.py`。

射击结果与参数：`docs/environment_reproduction/fixed_spin_kerr_cloud_shooting_20260917.json`。

独立 Leaver 对照：`docs/environment_reproduction/fixed_spin_cloud_frequency_20260917.json`。

本轮没有修改任何生产模块。
