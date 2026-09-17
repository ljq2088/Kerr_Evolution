# 标量变分与能流的逐式复核

本文件针对用户提供的 Dyson 等 PDF：第 3 页式 (13)–(19)，第 9–10 页式 (24)–(38)，并对照 Li 后续论文。这里区分恒等式、必须依赖背景方程的化简和论文印刷约定问题；不能把一处印刷问题自动认定为数值图偏差来源。

## 1. 从非密度 Klein–Gordon 算子重新变分

取度规号差 (−,+,+,+)，g_ab(λ)=g_ab+λh_ab，并令 h=g^{ab}h_ab、bar h_ab=h_ab−g_ab h/2。固定物理标量质量 μ。

\[
\delta g^{ab}=-h^{ab},\qquad
\delta\sqrt{-g}=\tfrac12h\sqrt{-g},\qquad
\delta\Gamma^c_{ab}=\tfrac12g^{cd}(\nabla_a h_{bd}+\nabla_bh_{ad}-\nabla_dh_{ab}).
\]

收缩连接变分：

\[
g^{ab}\delta\Gamma^c_{ab}=\nabla_a h^{ac}-\tfrac12\nabla^ch
=\nabla_a\bar h^{ac}\equiv C^c.
\]

对论文式 (14) 的非密度算子 Q=(□−μ²)Φ，有精确的一阶恒等式

\[
\boxed{\delta_g Q[\Phi]=-h^{ab}\nabla_a\nabla_b\Phi-C^b\nabla_b\Phi.}
\]

标量本身也变化时，δQ=(□−μ²)δΦ+δ_gQ。故一阶方程是

\[
(\Box-\mu^2)\delta\Phi=h^{ab}\nabla_a\nabla_b\Phi_c+C^b\nabla_b\Phi_c.
\]

洛伦兹规范 C^b=0 后得到 Dyson 式 (19)，与本地 `ThresholdCloud.lorenz_source` 的未迹反转 h 收缩相同。若错误地拿 bar h 收缩，则在精确背景上少了 μ²hΦc/2。迹反转只进入规范条件；不能因此把源张量也换成 bar h。

Dyson 印刷式 (17) 还包含 −hQ[Φc]/2。按其式 (14) 对 Q 的定义，此项不属于上述非密度算子变分。一个简单反例是常数 Φ、常数 μ：Q=−μ²Φ 与度规无关，精确 δ_gQ=0，但该附加项一般非零。若变分密度 sqrt(−g)Q，出现的是 +hQ/2，仍需注明算子定义改变。精确 Kerr 云 Q[Φc]=0 时，这一印刷差别消失，不能用它解释本轮 Kerr 图的显著差异。冻结复谱的虚部后，背景不再 on shell，必须重新说明采用哪一个近似方程，不能机械借用 on-shell 化简。

## 2. 轨道上的分布项

即使两侧真空区 C^b=0，拼接度规仍可能有

\[
C^b=[\bar h^{rb}]\delta(r-r_p).
\]

因此检查的是实际法向跳跃与云梯度的收缩，而非只看某一个度规分量：

\[
J_\delta^{\ell m}=2\pi\int_{-1}^1 S_{\ell m}^*\Sigma[\bar h^{rb}]\partial_b\Phi_c\,dx.
\]

本轮在两个新轨道 scalar22 上直接计算当前度规的单侧极限。L=6 时它对选定半径复场的相对改变量小于 1.6×10⁻⁹；L=8 控制仍小于 5.3×10⁻⁹。它不足以解释当前几十个百分点的场差。详见 `orbit_contact_scalar22_rp41p1.json` 与 `orbit_contact_scalar22_rp42p1.json`。这是有限 L 和所用外推的实际目标投影检验，不是完整度规匹配的无条件证明。

## 3. 分离变量和源权重

写 Φc=Rc(r)Sc(θ)exp(−iωct+imbφ)，单个度规傅里叶分量正比 exp(−imgΩpt+imgφ)，则

\[
m=m_b+m_g,\qquad\omega=\omega_c+m_g\Omega_p.
\]

复背景云不因取负 mg 而共轭；只有实度规满足 h_(−mg)=h_(mg)^*。把这一规则错用到整项 hΦc 会改变物理模式。

以 angular A 为分离常数：

\[
\frac1{\sin\theta}(\sin\theta S')'
+\left[a^2(\omega^2-\mu^2)\cos^2\theta-\frac{m^2}{\sin^2\theta}+A\right]S=0,
\]
\[
(\Delta R')'+\left[\frac{K^2}{\Delta}-\mu^2r^2-a^2\omega^2+2am\omega-A\right]R=J,
\quad K=(r^2+a^2)\omega-am.
\]

源 J 是 Σ 倍协变源的角向投影。本地名为 `lam` 的变量实际为 A；另一惯例 λ=A+a²ω²−2amω 必须先转换。计算是对 R 本身进行，若改用 u=√(r²+a²)R，径向算子、源、Wronskian 和边界振幅均须同时变换。

## 4. 用统一作用量系数核对归一化

写复标量作用量 L=−C(g^{ab}∂aΦ*∂bΦ+μ²|Φ|²)，则

\[
T_{ab}=C\left(\partial_a\Phi^*\partial_b\Phi+\partial_b\Phi^*\partial_a\Phi
-g_{ab}(|\nabla\Phi|^2+\mu^2|\Phi|^2)\right),
\quad j^a=-iC(\Phi^*\nabla^a\Phi-\Phi\nabla^a\Phi^*).
\]

在角函数单位球面范数和实频稳态假设下：

\[
N_\infty=2Ck|Z_\infty|^2,\quad
N_H=2C(r_+^2+a^2)(\omega-m\Omega_H)|Z_H|^2,
\]

其中传播支 k=sgn(ω)√(ω²−μ²)，束缚支无无穷远辐射通量。波的能量流为 ωN；云失去同样的 Noether 电荷时，轨道的有效损失为

\[
F_E^s=(\omega-\omega_c)N=m_g\Omega_pN,\qquad F_L^s=m_gN,
\qquad F_E^s=\Omega_pF_L^s.
\]

这解释了为什么代码必须区分 wave_energy 与 orbital_energy。Kerr 恒等式 r_+²+a²=2Mr_+ 也解释了视界面积因子；不能把 r_+² 直接当它使用。

Dyson 印刷 T_ab 式 (25) 与后面的电流/通量式有整体系数约定不一致；其 L 的印刷符号也应与显示的 KG 方程核对。必须把云质量归一化和通量中的 C 一并处理：固定 Mc 时，Φ∝C^(−1/2)，最终物理通量 C|Φ|² 不受纯约定的 C 影响。不能仅改通量或给场图乘一个拟合常数。

## 5. 不能忽略的原文数值说明

Dyson 第 4 页、式 (21) 后已经明确说明阈值前的小凹陷属于数值伪影，因为长波需要在更远处提取。这一段与本轮独立的长程库仑边界检验相符，因此该凹陷不能作为必须拟合的物理目标。但是这不自动证明 Li 的全部差异、或束缚场节点差异，都来自同一个数值设置。

本轮已分别完成真实源采样收敛、固定自旋准束缚背景和论文逐通道 tetrad 源检查。结果及适用范围见同目录总报告；它们仍不足以宣布所有场差已修复。


## 6. 径向 Green 核的测度与符号

令径向方程为 $(\Delta R')'+VR=J$，取 $u=R_{\rm in}$、$v=R_{\rm up}$，则

\[
W=\Delta(uv'-u'v),\qquad W'=0.
\]

使 $\Delta\partial_rG$ 在源点跳跃为 +1 的 Green 核是

\[
G(r,r')=\frac{u(r_<)v(r_>)}{W},\qquad
R(r)=\frac{v(r)}W\int_{r_H}^r u(r')J(r')\,dr'
+\frac{u(r)}W\int_r^{r_{\rm out}}v(r')J(r')\,dr'.
\]

这里 $J$ 已包含 $\Sigma$ 的角向投影，积分不再多乘 $\Delta$ 或 $\Sigma$。若径向方程预先除以 $\Delta$，其普通 Wronskian 也必须同时改为 $uv'-u'v=W/\Delta$，所得积分相同。只改其中一个约定会引入半径依赖的系统误差。

复合场需要先相干求和 $\sum R_{\ell m}S_{\ell m}e^{im\phi}$，再取模；不能把各通道的模相加。论文所画的 $\Phi^{(1,1)}$ 对应径向 $R$，不是中间变换变量 $\sqrt{r^2+a^2}R$。

## 7. 质量阈值附近，控制量不是单独的 |ω|

在 $M=1$ 时，令 $k^2=\omega^2-\mu^2$、$\beta=2\omega^2-\mu^2$。远区径向方程给出

\[
R\sim e^{ikr}r^{-1+i\beta/k}
\sim \frac{e^{ikr_*}r^{i\mu^2/k}}{\sqrt{r^2+a^2}}.
\]

第二种形式正是 Li 的边界定义；两种表示只有在 tortoise 的对数相位一起转换时才等价。传播支选择出射波；束缚支取 $k=i\kappa$、$\kappa>0$，变成 $e^{-\kappa r}r^{\beta/\kappa-1}$。普通有限阶 $1/r$ 展开的系数含 $1/k$ 的高次幂，因此在 $\omega\to\mu$ 时不均匀收敛。即使 $\omega\simeq0.3$ 并不很小，$|k|$ 也可能很小。

用 $y=\sqrt\Delta R$ 消去一阶导数，有

\[
y''+Q(r)y=0,\quad
Q=\frac{K^2+1-a^2}{\Delta^2}
-\frac{\mu^2r^2+A+a^2\omega^2-2am\omega}{\Delta}
=k^2+\frac{2\beta}{r}+O(r^{-2}).
\]

束缚模式的远区转折尺度约为 $2\beta/\kappa^2$。本例 $r_p=42.1$、标量22的精确外转折点约 $5169.706M$；在 $4000M$ 设置“已经处于衰减远区”的四阶条件会落在振荡区。此前用公开Notebook的精确四阶系数重算，确实发现该边界会改变场图，却不满足外边界收敛；$32000M\to64000M$ 的库仑边界控制则已稳定。

所以，MST/GSN 的选择不能单凭“ω小/大”自动解决这个 massive scalar 问题。必须另外控制 $k$、库仑参数 $\beta/k$、实际转折点与边界截断，检验同一物理解的 Wronskian、守恒通量和外边界收敛。这里保留通过收敛的库仑远区算法，不把失收敛设置当作贴图修复。

## 8. 印刷频率精度的独立敏感性控制

Li 显示 $\omega_c\simeq0.296294$，而独立复谱为 $0.29629353472561146+2.2166094\times10^{-9}i$。不能把显示的近似数字自动认定为作者计算输入。

仅固定旧源、更改 $a=0.88$ Green 频率的对照显示：在 $r_p=42.1$，使用该六位小数使标量22的两个径向振幅谷从 $54.3222,144.4219$ 移到 $59.3705,150.8266$；在 $r_p=41.1$，所检半径的最大复场改变量约 $6.59\times10^{-4}$。这体现了近束缚响应对频率的敏感性；不能根据零点附近的点态相对误差宣称全场都改变同样百分比。

全18模态的固定源六位数控制仍未匹配 Li 的图：$r_p=42.1$ 在 $r=50,100,150$ 的圆环幅度 RMS 偏差仍约 $48.5\%,13.2\%,44.7\%$。没有对频率进行参数寻优，也没有把较接近的结果当成正确答案。相应记录是 `scalar22_frequency_rounding_control.json`、`full_field_printed_frequency_fixed_source.json`；它们明确不是一致的新背景求解。

## 9. 对推导本身的独立数值检验

- 直接构造 $g\pm\epsilon h$、重算连接和非密度 KG 算子，再作有限差分，对照第1节变分恒等式；四个点包括 Schwarzschild/Kerr，最佳相对差 $1.12\times10^{-11}$ 至 $8.70\times10^{-11}$。测试使用 off-shell 标量和非 Lorenz 扰动，因此不会靠提前令背景方程或规范项为零掩盖错误。
- 从实际 L=6、mg=1 重构度规的一、二阶坐标导数直接组装 $\delta R_{ab}$，在 $r=3,20,50$ 离轨道点满足真空方程；按各连接变分项范数之和归一的残差为 $7.98\times10^{-16},1.05\times10^{-13},3.60\times10^{-12}$。$r=20$ 用Taylor阶数8/10对照稳定。纯规范 $h=\mathcal L_\xi g$ 的独立几何测试残差 $1.37\times10^{-16}$。
- 上述Ricci测试不直接使用重构所依据的 Teukolsky/GHP 方程作为检验算子，并已由另一审查者逐指标核对。但它仍是**离轨道局部真空测试**，不证明粒子 delta 源的全局幅度、完整轨道拼接或残余规范已与作者相同。

详见 `kg_metric_variation_audit.json` 和 `metric_vacuum_ricci_audit.json`。新增准束缚云及既有云/源/径向接口共24项测试通过。源公式遗漏及反事实的完整推导另见 `PAPER_DERIVATION_CONVENTION_AUDIT.md`；准束缚云的坐标、复角函数投影与归一化推导见 `GENERAL_KERR_QUASIBOUND_CLOUD.md`。
