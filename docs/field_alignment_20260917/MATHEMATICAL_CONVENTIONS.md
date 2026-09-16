# 本轮对齐使用的数学约定

以 M=1、α=μM=.3 为单位。背景是精确同步的 |211> 云，a=.8771530275949366，ωc=.29629324847975713。Li 等后续论文标注 a=.88、ωc≈.296294；这个参数差异在结果中保留，不能称严格相同背景。

归一化为 ζ=α³√(Mc/M)，物理扰动 δΦ=qζΦ11。求解文件把云质量归一到 Mc=1，因此场图使用缓存复场除以 α³。这不等于允许给源图或场图另拟合幅度常数；论文源图的绘图归一化应单独核验。

角函数满足

\[
\frac{1}{\sin\theta}\partial_\theta(\sin\theta\partial_\theta S)
+\left[a^2(\omega^2-\mu^2)\cos^2\theta-\frac{m^2}{\sin^2\theta}+A_{\ell m}\right]S=0,
\qquad 2\pi\int_0^\pi |S|^2\sin\theta\,d\theta=1.
\]

令 K=(r²+a²)ω−am、Δ=r²−2r+a²，则实际使用的径向方程为

\[
\partial_r(\Delta\partial_rR)+\left[\frac{K^2}{\Delta}-\mu^2r^2-a^2\omega^2+2am\omega-A_{\ell m}\right]R=J_{\ell m}.
\]

因此写成另一种常见形式时，λ=A+a²ω²−2amω，不能把已经含该移位的 λ 再代入代码的 A 位置。当前 `environment_radial.RadialGreen.lam` 的变量名虽然叫 lam，其存储含义是上式的 **A**。

Lorenz 条件消去一阶梯度项后，驱动源为

\[
J_{\ell m}(r)=2\pi\int_0^\pi S^*_{\ell m}(\theta)\,
\Sigma h^{ab}_{m_g}\nabla_a\nabla_b\Phi_c\,\sin\theta\,d\theta,
\qquad m=m_g+1,\quad\omega=\omega_c+m_g\Omega_p.
\]

这与 Li 径向源方程的 Σ 投影定义相同，不额外除 Δ 或 √(r²+a²)。论文源图若存在额外显示归一化，需要作者代码或独立明确约定，不能从一个拟合常数反推。

定义 W0=Δ(Rin Rup′−Rup Rin′)，连续响应是

\[
R(r)=R_{\rm up}(r)\int_{r_H}^{r}\frac{R_{\rm in}J}{W_0}dr'
+R_{\rm in}(r)\int_r^{\infty}\frac{R_{\rm up}J}{W_0}dr'.
\]

本轮另用独立插值和分裂 Gauss 积分复核这个表达式，没有通过改积分符号、负 m 相位或幅度来贴图。

场图是相干和的绝对值：

\[
\left|\Phi_{11}(r,\theta,\varphi,t)\right|
=\alpha^{-3}\left|\sum_{\ell,m}R_{\ell m}(r)S_{\ell m}(\theta)e^{i(m\varphi-\omega_m t)}\right|.
\]

本轮同用赤道 θ=π/2、t=0、X=r cosφ、Y=r sinφ。Li 图注将赤道写成 θ=0，与其角方程的极角定义矛盾；采用赤道的几何含义。图注“ell2..5全部模式”和正文“贡献无穷远通量的模式”也有歧义，因此分别输出：18个反射对称允许模式；6个正m模式 (2,2),(3,3),(4,2),(4,4),(5,3),(5,5)。后者在阈值之外仍保留m=2，符合正文所讨论的辐射到束缚转变。

本地视界轨道能量通量使用

\[
\dot E^{s,H}_{\ell m}=2(\omega-\omega_c)(r_+^2+a^2)(\omega-m\Omega_H)|Z_H|^2.
\]

后续论文明确指出原 Dyson 结果的分离常数和视界通量归一化有修正：[Li 等，v2，脚注1](https://arxiv.org/html/2507.02045v2)。这支持区分参考版本，不支持把每一项空间场差异都归咎于原论文。


## 近阈值外边界的独立控制

令 y=√Δ R，则长程齐次方程为

\[
y''+Q(r)y=0,\qquad
Q=k^2+\frac{2\beta}{r}+\frac{C_2}{r^2}+O(r^{-3}),
\quad \beta=2\omega^2-\mu^2,\quad C_2=12\omega^2-4\mu^2-a^2k^2-A.
\]

在束缚支 k=iκ，长程外转向点近似为

\[
r_t=\frac{\beta+\sqrt{\beta^2+\kappa^2C_2}}{\kappa^2}.
\]

本轮 rp=42.1 的 (2,2) 模给出 rt≈5169.72M。因此 4000M 仍可能处于振荡区，不能仅因其“远大于轨道半径”就使用短渐近截断。Whittaker 边界将 1/r 库仑尾重新求和；仍需外半径收敛，以控制遗漏的高阶势项。

零阶试验使用 R=e^(ikr) r^(−1+iβ/k)（忽略常相位）。在传播支，若误把有限半径系数当成单位无穷远振幅，真实流与名义流的比值为

\[
\frac{j_{\mathrm{actual}}}{j_{\mathrm{nominal}}}
=\frac{\Delta(R_b)}{R_b^2}
\left(1+\frac{\beta}{k^2R_b}\right).
\]

rp=41.1、Rb=4000M 时给出 1.493899，与直接求守恒流相同。这种通量偏差不是常规 ODE 容差能够修复的。

Li 正文采用的前因子实际上是 exp(ikr*) r^(iμ²/k)/√(r²+a²)，并保留四阶 1/r 展开；这和上述零阶实验的有限半径数据不能混同。公开仓库只给出边界系数，没有给出原图的实际外截断半径：[作者 Boundary Coefficients.nb](https://github.com/dongjun826/EMRI-in-scalar-clouds)。本实验不能确定作者的真实积分设置。


按 Li 的精确前因子 P，令 R=P∑(B_j/r^j)，使用

\[
\frac{P'}P=ik\frac{r^2+a^2}{\Delta}+\frac{i\mu^2}{kr}-\frac{r}{r^2+a^2}
\]

代入径向方程独立递推 B1…B4，与公开 notebook 的系数在两个算例中以 80 位精度比对，最大相对差 1.43×10⁻⁴⁵。四阶展开置于 4000M 时，束缚场的 r=50、150 圆周分布比库仑结果更接近参考图，但移动至 8000M、32000M 后又回到库仑结果。这验证了边界截断敏感性；没有证明原文实际采用的截断位置。
