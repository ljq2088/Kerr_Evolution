# A6 连续演化方程与数值接口

- 主文：[DERIVATION.md A6 section](../DERIVATION.md)
- 连续 operator：[附录 D.4](W02_1_4_FIXED_M_Y_DERIVATION.md#a4-y-operator)
- 数值 source：[A5 四块公式](a5_four_block_source_terms.md)
- 主要依据：`P13` Eqs. (22)--(25)、Sec. V 与 Appendix C；axis regularity 参考 `P18` 并由本文 indicial equation 独立推导
- 状态：A6 continuous contract 与三组 targeted checks 已获 owner 接受；boundary-adjacent finite-difference closure 未接受；未作 full-endpoint claim
- 范围：连续 PDE、一阶化、worldline/stage interface、continuum characteristics、axis/SCRI projection contract
- 非结论：transition prescription、finite-difference matrices、CFL factor、startup 和 production admission 仍是 W02.3 工程工作

## Post-ISCO 初值接口

### ISCO quantities

对 baseline $0\le\chi=a/M<1$ 的 prograde equatorial orbit，令 $\tilde r=r/M$，

$$
Z_1=1+(1-\chi^2)^{1/3}
\left[(1+\chi)^{1/3}+(1-\chi)^{1/3}\right],
$$

$$
Z_2=\sqrt{3\chi^2+Z_1^2},
$$

$$
\boxed{
\tilde r_{\rm ISCO}
=3+Z_2-\sqrt{(3-Z_1)(3+Z_1+2Z_2)}
}.
$$

取 $v=\tilde r_{\rm ISCO}^{-1/2}$，

$$
\boxed{
E_{\rm ISCO}
=\frac{1-2v^2+\chi v^3}
{\sqrt{1-3v^2+2\chi v^3}}
},
$$

$$
\boxed{
\frac{L_{z,\rm ISCO}}M
=\frac{1-2\chi v^3+\chi^2v^4}
{v\sqrt{1-3v^2+2\chi v^3}},
\qquad
M\Omega_{\rm ISCO}
=\frac1{\tilde r_{\rm ISCO}^{3/2}+\chi}
}.
$$

精确 ISCO circular geodesic 不会自行 plunge。A6 的实际 continuous interface 因此接收

$$
\boxed{
r_0=r_{\rm ISCO}-\delta r,
\qquad
E_{\rm p}=E_{\rm ISCO}-\delta E,
\qquad
L_{z,\rm p}=L_{z,\rm ISCO}-\delta L_z
},
$$

以及产生 $(\delta r,\delta E,\delta L_z)$ 的 transition manifest。必须检查

$$
\boxed{
r_+<r_0<r_{\rm ISCO},
\qquad
\mathcal R(r_0;E_{\rm p},L_{z,\rm p})\ge0,
\qquad
E_{\rm p}-\Omega_HL_{z,\rm p}>0
}.
$$

Ori--Thorne、flux truncation、$q_{\rm tr}$ 与 match prescription 是待 owner/W02.3 决定的 provenance，不由以上代数唯一确定。

## 唯一 $T$-parameterized worldline

只积分

$$
\boxed{
z_p(T)=\bigl(R_p(T),0,\Phi_p(T)\bigr)
},
$$

$$
\boxed{
\dot R_p=\frac{u_p^R}{u_p^T},
\qquad
\dot\Phi_p=\frac{u_p^\Phi}{u_p^T}
},
$$

并取

$$
T_0=0,
\qquad
R_p(0)=\frac{L^2}{r_0},
\qquad
\Phi_p(0)=0.
$$

每次 RHS evaluation 用

$$
r_p(T)=\frac{L^2}{R_p(T)}
$$

评价 BL radial functions 和 regular hyperboloidal velocities。恒等式

$$
\frac{dr_p}{dT}
=-\frac{r_p^2}{L^2}\dot R_p
=-\frac{r_p^2u_p^R}{L^2u_p^T}
$$

是 diagnostic，不定义第二条轨迹。BL $t_p,\phi_p$ 不积分；如需诊断，只按

$$
t_p=T-h(r_p),
\qquad
\phi_p=\Phi_p-q(r_p)
$$

重建。Trajectory contract 输出 dense stage-time interpolant、A5 所需 jets、horizon event 和 Gaussian exterior tail 的 source-exit event。

<a id="a6-pde"></a>

## 完整二阶 PDE：A--F

记

$$
d(y)=1-y^2,
$$

$$
\mathcal A_{-2,m}[\psi]
=d\psi_{,yy}-2y\psi_{,y}
-\frac{(m+2y)^2}{d}\psi-2\psi.
$$

完整 equation 写成

$$
\boxed{
A\psi_{,TT}+B\psi_{,TR}+C\psi_{,RR}
-\mathcal A_{-2,m}[\psi]
+D\psi_{,T}+E_R\psi_{,R}+F\psi
=S_m^{\rm num}
},
$$

其中 $\psi=\psi_{4,m}$，

$$
\boxed{
A=8M\left(2M-\frac{a^2R}{L^2}\right)
\left(1+\frac{2MR}{L^2}\right)-a^2d
},
$$

$$
\boxed{
B=-2\left[
L^2-(8M^2-a^2)\frac{R^2}{L^2}
+\frac{4a^2MR^3}{L^4}
\right]
},
$$

$$
\boxed{
C=-\left(L^2-2MR+\frac{a^2R^2}{L^2}\right)
\frac{R^2}{L^2}
},
$$

$$
\boxed{
\begin{aligned}
D={}&2iam\left(1+\frac{4MR}{L^2}\right)\\
&+2\left[
2M\left(2-3\frac{a^2R^2}{L^4}\right)
-\frac{a^2R}{L^2}+2iay
\right],
\end{aligned}
},
$$

$$
\boxed{
E_R=\frac{2iamR^2}{L^2}
+2R\left[
1+\frac{MR}{L^2}-2\frac{a^2R^2}{L^4}
\right]
},
$$

$$
\boxed{
F=\frac{2iamR}{L^2}
-2\frac{MR}{L^2}-2\frac{a^2R^2}{L^4}
}.
$$

这些 coefficients 是附录 D.4 取 $s=-2$ 后的逐项 regrouping；不是新的方程选择。

## 连续 P/Q 一阶系统

按 `P13` Eq. (23) 定义

$$
\boxed{
P=A\psi_{,T}+B\psi_{,R}+D\psi,
\qquad
Q=\psi_{,R}
}.
$$

由于 coefficients 与 $T$ 无关，二阶 PDE 等价于

$$
\boxed{
\psi_{,T}=\frac{P-BQ-D\psi}{A}
},
$$

$$
\boxed{
P_{,T}
=S_m^{\rm num}-C Q_{,R}
+\mathcal A_{-2,m}[\psi]
-E_RQ-F\psi
},
$$

$$
\boxed{
Q_{,T}=\partial_R
\left(\frac{P-BQ-D\psi}{A}\right)
}.
$$

Constraint 是

$$
\mathcal C_Q=Q-\psi_{,R}=0.
$$

连续三场系统是理论恒等式。Method-of-lines baseline 则采用 constrained two-field state $(P,\psi)$，每个 stage 重算

$$
Q=D_R^{(1)}\psi,
\qquad
Q_{,R}=D_R^{(2)}\psi.
$$

这里 $D_R^{(2)}$ 是独立 second-derivative matrix，不能预设为 $(D_R^{(1)})^2$。Source 只以 $+S_m^{\rm num}$ 进入 $P_{,T}$。

## 从 field initial data 到 RK4 stage

Initial field interface 是

$$
\psi_0(R,y)=\psi(T_0,R,y),
\qquad
\Pi_0(R,y)=\partial_T\psi(T_0,R,y),
$$

$$
\boxed{
P_0=A\Pi_0+B D_R^{(1)}\psi_0+D\psi_0
}.
$$

Zero-field $(\psi_0,\Pi_0)=(0,0)$ 只是一项 startup candidate。动态数组为

$$
U(T)=\begin{pmatrix}P\\ \psi\end{pmatrix}.
$$

任意 stage $(T_s,P_s,\psi_s)$ 的数据流是

$$
Q_s=D_R^{(1)}\psi_s,
\qquad
Q_{R,s}=D_R^{(2)}\psi_s,
\qquad
\mathcal A_{s}^{\rm FD}
=\mathcal A_{-2,m}^{\rm FD}[\psi_s],
$$

$$
z_p(T_s)
\longrightarrow
\text{worldline/projection/Gaussian jets}
\longrightarrow
\{S_{m,A}^{\rm num}(T_s)\}_{A=1}^4
\longrightarrow S_m^{\rm num}(T_s),
$$

$$
\boxed{
\dot\psi_s=\frac{P_s-BQ_s-D\psi_s}{A}
},
$$

$$
\boxed{
\dot P_s=S_m^{\rm num}(T_s)-C Q_{R,s}
+\mathcal A_s^{\rm FD}-E_RQ_s-F\psi_s
}.
$$

记该 RHS 为 $\dot U=\mathcal F(T,U)$。Classical RK4 使用

$$
k_1=\mathcal F(T_n,U_n),
$$

$$
k_2=\mathcal F\left(T_n+\frac{\Delta T}{2},
U_n+\frac{\Delta T}{2}k_1\right),
$$

$$
k_3=\mathcal F\left(T_n+\frac{\Delta T}{2},
U_n+\frac{\Delta T}{2}k_2\right),
$$

$$
k_4=\mathcal F(T_n+\Delta T,U_n+\Delta T k_3),
$$

$$
\boxed{
U_{n+1}=U_n+\frac{\Delta T}{6}
(k_1+2k_2+2k_3+k_4)
}.
$$

每次 $\mathcal F$ 调用都必须在对应 stage time 重采样 trajectory 并重算 source；不能只在完整 time-step 起点更新粒子。

## Characteristics 与 radial boundaries

Radial principal polynomial 是

$$
A v_R^2-Bv_R+C=0,
$$

故

$$
\boxed{
v_R^\pm
=\frac{B/2\pm\sqrt{B^2/4-AC}}{A}
}.
$$

Angular principal speed magnitude 是

$$
\boxed{
|v_y|=\sqrt{\frac{1-y^2}{A}}
}.
$$

在 SCRI+，$R=0$，

$$
C=0,
\qquad B=-2L^2,
$$

$$
\boxed{
(v_R^+,v_R^-)_{\mathscr I^+}
=\left(0,-\frac{2L^2}{A(0,y)}\right)
}.
$$

一个 characteristic tangent，另一个向 decreasing $R$ 离开计算域；没有从 $R<0$ 指向 interior 的 incoming field characteristic。

在 future horizon $R=R_H=L^2/r_+$，

$$
C=0,
\qquad B=4MR_H>0,
$$

$$
\boxed{
(v_R^+,v_R^-)_{\mathcal H^+}
=\left(\frac{4MR_H}{A(R_H,y)},0\right)
}.
$$

一个 characteristic tangent，另一个向 increasing $R$ 穿出 exterior domain；同样没有 incoming physical field data。One-sided radial rows 是数值 closure，不是额外物理 boundary condition。

## Axis indicial equation 与 local factor-aware closure

对 $s=-2$，

$$
\mathcal A_{-2,m}[\psi]
=d\psi_{,yy}-2y\psi_{,y}
-\frac{(m+2y)^2}{d}\psi-2\psi,
\qquad d=1-y^2.
$$

表面的发散来自最后一项中的 $d^{-1}$。关键不是单独删除这一项，而是先确定 smooth spin-weighted field 在轴上允许的幂次，再证明完整 angular operator 中相同阶的奇异项彼此抵消。

### North axis：$y=-1$

令

$$
x=1+y,\qquad y=x-1,
$$

则 north axis 对应 $x=0$，并且

$$
d=1-y^2=x(2-x),
\qquad
m+2y=(m-2)+2x.
$$

记

$$
a_-=m-2.
$$

在 north axis 附近作局部分解

$$
\boxed{
\psi(y)=x^p g(x)
},
$$

其中 $g(x)$ 在 $x=0$ 有普通 Taylor expansion，且不预先假定 $p$。因为 $\partial_y=\partial_x$，

$$
\psi_{,y}
=p x^{p-1}g+x^p g',
$$

$$
\psi_{,yy}
=p(p-1)x^{p-2}g+2p x^{p-1}g'
+x^p g''.
$$

先代入前两个微分项：

$$
\begin{aligned}
d\psi_{,yy}-2y\psi_{,y}
=x^p\Bigg[
&x(2-x)g''
+\left(2p(2-x)-2(x-1)\right)g'\\
&+\left(\frac{2p^2}{x}-p(p+1)\right)g
\Bigg].
\end{aligned}
$$

势项可精确拆成

$$
\frac{(a_-+2x)^2}{x(2-x)}
=\frac{a_-^2}{2x}
+\frac{a_-^2/2+4a_-+4x}{2-x}.
$$

所以完整算符为

$$
\begin{aligned}
\mathcal A_{-2,m}[x^p g]
=x^p\Bigg\{
&x(2-x)g''
+\left[2p(2-x)-2(x-1)\right]g'\\
&+\left[
\frac{2p^2-a_-^2/2}{x}
-p(p+1)-2
-\frac{a_-^2/2+4a_-+4x}{2-x}
\right]g
\Bigg\}.
\end{aligned}
$$

这里唯一的 negative power 是

$$
\frac{2p^2-a_-^2/2}{x}g.
$$

要求 smooth field 经过完整算符后不产生比 $\psi$ 更强的 $x^{-1}$ 发散，必须令

$$
2p^2-\frac{a_-^2}{2}=0.
$$

取 regular branch $p\ge0$，

$$
\boxed{
p_-=p=\frac{|a_-|}{2}
=\frac{|m-2|}{2}
}.
$$

代入后得到没有任何 $x^{-1}$ 项的 north regular operator：

$$
\boxed{
\begin{aligned}
\mathcal A_{-2,m}[x^{p_-}g]
=x^{p_-}\Bigg\{
&x(2-x)g''
+\left[2p_-(2-x)-2(x-1)\right]g'\\
&+\left[
-p_-(p_-+1)-2
-\frac{(m-2)^2/2+4(m-2)+4x}{2-x}
\right]g
\Bigg\}.
\end{aligned}
}
$$

大括号内的所有 coefficients 在 $x=0$ 都有限。特别地，

$$
\left.
x^{-p_-}\mathcal A_{-2,m}[x^{p_-}g]
\right|_{x=0}
=(4p_-+2)g'(0)
+\left[-2p_-^2-p_--2-2(m-2)\right]g(0).
$$

例如 $m=2$ 时 $p_-=0$，

$$
\boxed{
\left.\mathcal A_{-2,2}[\psi]\right|_{y=-1}
=2g'(0)-2g(0)
},
$$

所以 north axis 上的 field 可以有限且非零，并不存在真实发散。

### South axis：$y=+1$

令

$$
z=1-y,\qquad y=1-z,
$$

则 south axis 对应 $z=0$，

$$
d=z(2-z),
\qquad
m+2y=(m+2)-2z.
$$

记

$$
a_+=m+2,
$$

并写

$$
\boxed{
\psi(y)=z^p g(z)
}.
$$

注意此时

$$
\partial_y=-\partial_z,
\qquad
\partial_y^2=\partial_z^2.
$$

因此

$$
\psi_{,y}
=-p z^{p-1}g-z^p g',
$$

$$
\psi_{,yy}
=p(p-1)z^{p-2}g+2p z^{p-1}g'
+z^p g''.
$$

微分部分成为

$$
\begin{aligned}
d\psi_{,yy}-2y\psi_{,y}
=z^p\Bigg[
&z(2-z)g''
+\left(2p(2-z)+2(1-z)\right)g'\\
&+\left(\frac{2p^2}{z}-p(p+1)\right)g
\Bigg].
\end{aligned}
$$

而

$$
\frac{(a_+-2z)^2}{z(2-z)}
=\frac{a_+^2}{2z}
+\frac{a_+^2/2-4a_++4z}{2-z}.
$$

于是

$$
\begin{aligned}
\mathcal A_{-2,m}[z^p g]
=z^p\Bigg\{
&z(2-z)g''
+\left[2p(2-z)+2(1-z)\right]g'\\
&+\left[
\frac{2p^2-a_+^2/2}{z}
-p(p+1)-2
-\frac{a_+^2/2-4a_++4z}{2-z}
\right]g
\Bigg\}.
\end{aligned}
$$

唯一的 negative power 是

$$
\frac{2p^2-a_+^2/2}{z}g,
$$

所以 regular branch 是

$$
\boxed{
p_+=\frac{|a_+|}{2}
=\frac{|m+2|}{2}
}.
$$

代入后，

$$
\boxed{
\begin{aligned}
\mathcal A_{-2,m}[z^{p_+}g]
=z^{p_+}\Bigg\{
&z(2-z)g''
+\left[2p_+(2-z)+2(1-z)\right]g'\\
&+\left[
-p_+(p_++1)-2
-\frac{(m+2)^2/2-4(m+2)+4z}{2-z}
\right]g
\Bigg\},
\end{aligned}
}
$$

其大括号内在 $z=0$ 同样有限，并且

$$
\left.
z^{-p_+}\mathcal A_{-2,m}[z^{p_+}g]
\right|_{z=0}
=(4p_++2)g'(0)
+\left[-2p_+^2-p_+-2+2(m+2)\right]g(0).
$$

### $m=2,4$ 的结果

$$
\boxed{
\begin{array}{c|cc}
m&p_-(y=-1)&p_+(y=+1)\\ \hline
2&0&2\\
4&1&3
\end{array}
}.
$$

因此

$$
\psi_{4,2}\sim
\begin{cases}
g_-(1+y),&y\to-1,\\
(1-y)^2g_+(1-y),&y\to+1,
\end{cases}
$$

$$
\psi_{4,4}\sim
\begin{cases}
(1+y)g_-(1+y),&y\to-1,\\
(1-y)^3g_+(1-y),&y\to+1.
\end{cases}
$$

这里的 $g_\pm$ 都是 smooth functions。表面上的 $d^{-1}$ 没有被忽略；它的 $x^{-1}$ 或 $z^{-1}$ 部分恰好与 $d\psi_{,yy}-2y\psi_{,y}$ 产生的同阶项相消。

### “重定义”究竟做了什么

$$
\psi=f_\pm g_\pm,
\qquad
f_-=(1+y)^{p_-},
\qquad
f_+=(1-y)^{p_+}
$$

首先是一种 **局部解析分解**，用来显式提出 smooth spin-weighted field 已知的轴上零点；它不是删除原方程中的奇异系数。上面的计算证明

$$
\boxed{
\mathcal A_{-2,m}[f_\pm g_\pm]
=f_\pm\,\widetilde{\mathcal A}_\pm[g_\pm]
},
$$

其中 $\widetilde{\mathcal A}_\pm$ 的 coefficients 在相应 axis 上有限。换言之，

$$
\frac{\mathcal A_{-2,m}[f_\pm g_\pm]}{f_\pm}
$$

有有限的 axis limit。发散消失的原因是两条明确的 indicial conditions

$$
4p_-^2=(m-2)^2,
\qquad
4p_+^2=(m+2)^2
$$

精确消除了唯一的 negative Laurent coefficient，而不是因为数值上把 $1/d$ 截断为有限值。

当前 baseline 不在全域演化 $g_\pm$。它仍然演化原变量 $\psi$，只在靠近相应 axis 的 stencil 中使用 $\psi=f_\pm g_\pm$ 的局部结构。

### Open grid 与 local factor-aware derivative rows

Angular grid 排除 $y=\pm1$，因此 interior nodes 上 $d>0$，原算符本身可以直接评价。问题只是在靠近 axis 时，普通 polynomial stencil 未必保持上述大项之间的解析 cancellation。为此，对 north/south boundary-adjacent rows 分别使用

$$
f_-(y)=(1+y)^{p_-},
\qquad
f_+(y)=(1-y)^{p_+}.
$$

在 stencil center $y_i$ 附近，把 smooth quotient 近似为

$$
g_\pm(y)
\approx\sum_{k=0}^{N_{\rm st}-1}
c_k(y-y_i)^k,
$$

所以

$$
\psi(y)
\approx\sum_{k=0}^{N_{\rm st}-1}
c_k f_\pm(y)(y-y_i)^k.
$$

要求离散 $n$ 阶 derivative 对每个 basis function 都精确，得到 weights equations

$$
\boxed{
\sum_j w_{ij}^{(n,\pm)}
f_\pm(y_j)(y_j-y_i)^k
=\left.
\partial_y^n\left[
f_\pm(y)(y-y_i)^k
\right]\right|_{y=y_i},
\quad
k=0,\ldots,N_{\rm st}-1
}.
$$

解这个 $N_{\rm st}\times N_{\rm st}$ linear system 后，

$$
\partial_y^n\psi(y_i)
\approx\sum_jw_{ij}^{(n,\pm)}\psi(y_j).
$$

这与普通 finite difference 的区别仅在于 basis 从 $(y-y_i)^k$ 换成了 $f_\pm(y)(y-y_i)^k$。Interior rows 仍使用普通 finite differences，动态变量仍是 $\psi$。

Endpoint reconstruction 使用相同的 local polynomial：

$$
\boxed{
\psi(y_\pm)=0\quad(p_\pm>0)
}.
$$

若 $p_\pm=0$，则 $f_\pm=1$，所以

$$
\boxed{
\psi(y_\pm)=g_\pm(y_\pm)
}
$$

由邻近 open nodes 对 smooth $g_\pm$ 外推得到。对于当前 modes，只有 $m=2$ 的 north axis 属于 $p_-=0$；其余三个 mode/axis combinations 的 endpoint field 都严格为零。

实现时不得在 endpoint 分别计算

$$
d\psi_{,yy},\qquad
-2y\psi_{,y},\qquad
-\frac{(m+2y)^2}{d}\psi.
$$

若确实需要 endpoint angular-operator value，必须使用上面已经消去 negative power 的 $\widetilde{\mathcal A}_\pm[g_\pm]$ finite limit。实际 open nodes、stencil width 和 production endpoint rows 仍属于 W02.3，并必须用 analytic spin-weighted harmonics 检验。Global field redefinition 只保留为 local closure 失败时的 fallback。

## Gaussian spatial strategy

候选 baseline 是

$$
\boxed{
\sigma_R=4\Delta R,
\qquad
\sigma_y=4\Delta y
}.
$$

“四个 cells”在 refinement 中固定，因此 physical widths 随 grid spacing 缩小。这是 coupled grid--Gaussian convergence strategy，不是从连续方程推导出的唯一宽度，也不能单独给出 fixed-source stencil order。异常时再做 independent width scan。

## SCRI+ projection

在 $R=0$，$T=u$，连续 projection 是

$$
\boxed{
\psi_{4,\ell m}(T)
=2\pi\int_{-1}^{1}
{}_{-2}\!\bar Y_{\ell m}(\arccos(-y),0)
\psi_{4,m}(T,0,y)\,dy
}.
$$

若 angular nodes 是 $y_j$，quadrature weights 是 $w_j$，定义 projection row

$$
\boxed{
\mathsf P_j^{(\ell m)}
=2\pi w_j\,{}_{-2}\!\bar Y_{\ell m}
(\arccos(-y_j),0)
},
$$

$$
\boxed{
\psi_{4,\ell m}(T)
\approx\sum_j\mathsf P_j^{(\ell m)}
\psi_{4,m}(T,0,y_j)
}.
$$

$\mathsf P$ 是 quadrature/projection matrix，不是 derivative matrix。Open-node quadrature 或先重建 endpoints 的 quadrature 都可作为候选，但必须检查 analytic target-mode recovery、相邻 $\ell$ leakage 和 normalization。权威输出仍是 complex $\psi_4$；strain conversion 是另行受控的派生量。

## Theory contract 与工程开放项

### 已推导并可进入 A6 checks

- ISCO formulas 与 post-ISCO input inequalities；
- 唯一 $(R_p,\Phi_p)$ trajectory equations；
- 完整 A--F 二阶 PDE；
- P/Q continuous first-order equivalence；
- initial-field 到每个 RK4 stage 的 data flow；
- radial/angular principal speeds 与 continuum boundary directions；
- axis indicial powers、factor-aware consistency equation 与 endpoint rule；
- continuous/discrete SCRI+ projections。

### 仅为工程策略

- transition/flux/match prescription 和 $q_{\rm tr}$；
- zero-field 或其他 startup；
- grid sizes、stencil width、open-node placement 和 matrix construction；
- CFL factor、output cadence 和 source-inactive policy；
- four-cell Gaussian preset 的 production admission；
- coupled grid--Gaussian convergence 和 frequency-domain benchmark。

### 最小 A6 Mathematica checks

1. 从附录 D.4 独立抽取 $s=-2$ coefficients，与 $A$--$F$ 逐项比较，并把 $P,Q$ definitions 代回恢复完整二阶 PDE。
2. 对 $v_R^\pm$ 在 $R=0,R_H$ 取极限，核对 $C=0$、$B=-2L^2$ 与 $B=4MR_H$ 及 characteristic directions。
3. 对 north/south axis 重做 indicial residual，并对 $m=2,4$ 核对 $(p_-,p_+)=(0,2),(1,3)$；用小型 polynomial basis 检查 factor-aware weight consistency equation。
4. 用 analytic spin-weighted harmonics 检查 continuous/discrete projection normalization 和相邻 $\ell$ leakage。Axis finite-limit matrix 的具体 coefficients 依赖 W02.3 nodes/stencil，不在 A6 Mathematica 中虚构。

A6 targeted results 见 [summary](../../../../data/kerr_point_particle_verification/results/a6/a6_targeted_summary.json)：`FirstOrderIdentityChecked`、`CharacteristicsChecked`、`AxisProjectionChecked` 均为 `true`，且 `FullEndpointClaimMade=false`。Owner 已接受这一连续/符号范围，但 boundary-adjacent finite-difference rows 仍需在 W02.3 建立理解和实际数值证据；本文不声称数值 stability、convergence 或 production readiness。
