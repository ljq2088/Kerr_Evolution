# Kerr hyperboloidal minimal gauge：方程与代码对应

## 1. 文献与记号

以 Zotero `8AW4IT3C` 的综述第 3.2 节为框架入口，以同收藏夹 `R7UUGT4R`
的 Ripley (2022), arXiv:2202.03837v2，第 2 节、式 (5)–(11) 为 Kerr 实现依据。
另核对综述参考文献 [113]：Macedo (2020), arXiv:1910.13452v2，
第 4.1 节的 radial function fixing minimal gauge。

下文 `a` 是有量纲 Kerr 旋转参数，`chi=a/M` 是命令行的 `--spin`；
标量场的 Teukolsky 自旋权重为零，不要与 `a` 混淆。
`L=r_+=M+sqrt(M^2-a^2)`，`Delta=r^2-2Mr+a^2`，`Sigma=r^2+a^2 cos^2(theta)`。

## 2. 坐标与场

从 Boyer–Lindquist 坐标出发：

\[
dv=dt+\frac{r^2+a^2}{\Delta}dr,\quad
d\varphi=d\phi+\frac a\Delta dr,\quad
\tau=v-2r-4M\ln(r/L),\quad \sigma=L/r.
\]

Ripley 的 ingoing 时间是 `v_R=v-r`，因此其高度函数导数为 `-1-4M/r`，
与这里的 `tau=v_R-r-4M ln(r/L)` 完全一致（允许常数时间平移）。
Ripley 使用 `rho=1/r`，这里用 `sigma=L rho` 把视界固定在 1。
Macedo 使用无量纲时间，本文代码保留有量纲时间，关系为 `tau_Macedo=tau/L`。

方位角保持 ingoing `varphi`，**不能**继续用未变换的 Boyer–Lindquist 方位角套用系数。
在 `a=0` 极限，`tau=t+r_*-2r-4M ln(r/L)`，恢复综述第 3.2 节的 minimal gauge。
在 `a!=0` 时高度函数形式仍简单，但方程系数与角向耦合必须使用完整 Kerr 度规。

设 `u=r Phi`，有

\[
u(\tau,\sigma,\theta,\varphi)=\sum_{\ell=|m|}^{\ell_{\max}}
u_{\ell m}(\tau,\sigma)Y_{\ell m}(\theta,\varphi).
\]

此变换是共形场正则化；在零无穷，`Phi` 本身趋零，有限的 `u` 才是辐射场。
球谐采用 `exp(im varphi)` 约定。固定 `m` 扇区独立，固定 `l` 在 Kerr 上并不独立。

## 3. 独立链式推导

ingoing `(v,r,theta,varphi)` 中的标量方程乘以 `Sigma` 为

\[
a^2\sin^2\theta\Phi_{vv}+2(r^2+a^2)\Phi_{vr}
+2a\Phi_{v\varphi}+\Delta\Phi_{rr}+2a\Phi_{r\varphi}
+2r\Phi_v+\Delta'\Phi_r+\Delta_{S^2}\Phi=0.
\]

这里 `Delta_{S^2}` 是完整单位球面 Laplacian，包括方位角项。
代入 `Phi=u/r`，再用
`partial_r|v = -(sigma^2/L) partial_sigma - (2+4M sigma/L) partial_tau`。
化简后各系数为多项式，无需在端点对无穷大数相消。

定义 `H=2+4M sigma/L` 及单位矩阵 `I`，得到矩阵方程

\[
A\ddot{\boldsymbol u}+C\partial_\sigma\dot{\boldsymbol u}
+D\partial_\sigma^2\boldsymbol u+E\dot{\boldsymbol u}
+F\partial_\sigma\boldsymbol u+V\boldsymbol u=0,
\]

\[
\begin{aligned}
A&=H(-8M^2+4Ma^2\sigma/L)I+a^2 S,\\
C&=2L+(2a^2-16M^2)\sigma^2/L+8Ma^2\sigma^3/L^2,\\
D&=\sigma^2(1-2M\sigma/L+a^2\sigma^2/L^2),\\
E&=(2a^2-16M^2)\sigma/L+12Ma^2\sigma^2/L^2+2iam(1-H),\\
F&=2\sigma-6M\sigma^2/L+4a^2\sigma^3/L^2-2iam\sigma^2/L,\\
V_{\ell\ell'}&=\left[-2M\sigma/L+2a^2\sigma^2/L^2-2iam\sigma/L
-\ell(\ell+1)\right]\delta_{\ell\ell'}.
\end{aligned}
\]

`C,D,E,F` 对所有 `l` 相同；`S_ll'=<Y_lm,sin^2(theta)Y_l'm>`。
以上为从 Klein–Gordon 方程的独立推导；将整体乘 `-1` 并令 `sigma=L rho`，
与 Ripley 式 (11) 在场自旋权重零时逐项一致，测试对多组质量和旋转参数验证了这一映射。

## 4. 角向与端点

用 `cos(theta)Y_lm = c_(l+1)Y_(l+1,m)+c_l Y_(l-1,m)`，
`c_l=sqrt((l^2-m^2)/(4l^2-1))` 构造 `S=I-X^2` 的投影。
必须先包含中间模态 `lmax+1` 再取子矩阵，否则最顶端的对角元素错误。
本实现对该投影另用 Gauss–Legendre 积分检查。

`A` 在计算域负定，因此 `tau=constant` 切片是类空的。
`D=sigma^2(1-sigma)(1-(a/L)^2 sigma)` 在视界及零无穷为零。
对径向主部 `A u_tt+C u_ts+D u_ss`，特征速度满足
`A c^2-C c+D=0`。端点一支速度为零，另一支为 `C/A`：
在零无穷为负（流出 `sigma>=0`），在视界为正（流出 `sigma<=1`）。
因此在端点直接使用退化后的演化方程，不施加 Dirichlet、Neumann 或 Sommerfeld 外加数据。
代码域终止于事件视界；切片可以穿越视界，但本次没有演化黑洞内部。

## 5. 时间演化与精度

径向采用 Chebyshev–Lobatto 配点及一、二阶微分矩阵。
设 `p=u_tau`，每一径向点按
`u_tau=p`, `p_tau=-A^{-1}(C D_sigma p+D D_sigma^2 u+E p+F D_sigma u+V u)`
组成常微分方程组，由 SciPy DOP853 自适应积分。
没有添加数值耗散或人工吸收层。显式方法仍受空间离散稳定性约束，高分辨率时更慢；
收敛检查必须同时改变径向点数、角向截断和时间容差。

初始 `u` 在单个球谐上取
`[sigma(1-sigma)/(c(1-c))]^2 exp(-((sigma-c)/w)^2)`，`p=0`。
该函数在计算闭区间解析，两端函数及一阶导数为零；不是严格紧支撑于开区间的 bump。
默认 `c=.45,w=.12`。本次结果不能直接用于声称已经测准 Price 尾指数或 QNM 频率。
