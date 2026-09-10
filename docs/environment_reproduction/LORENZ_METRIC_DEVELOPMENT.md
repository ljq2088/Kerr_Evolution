# 非静态 Lorenz 度规：推导与验证边界

这是开发记录，不是已验收的点粒子度规。采用 M=1、Fourier 因子
exp(-i omega t+i m phi)，仅实现 omega 非零、r 不等于轨道半径的区域。
参考为 [2406.12510v3](https://arxiv.org/html/2406.12510v3)。

## 从曲率到张量

论文的复 AAB 场在真空区域为

\[
h^{\rm AAB}_{ab}=\frac43\left[
 (\mathcal S_4^\dagger\zeta^4\psi_0)_{ab}
 -(\mathcal S_0^\dagger\zeta^4\psi_4)_{ab}\right].
\]

令 DKW 规范向量为 \(\xi^{\rm DKW}\)，则自旋2部分为

\[
h^{(2)}_{ab}=\frac{i}{\omega}
 \left[h^{\rm AAB}_{ab}-2\nabla_{(a}\xi^{\rm DKW}_{b)}\right].
\]

这里使用原始曲率幅度，不能把公开表的时间积分幅度再代入同一表达式。
二形式散度的指标方向在代码中显式记录；改变方向会破坏 Lorenz 条件。

独立验证从 Christoffel 符号开始，而不是再次调用重构算符：

\[
\delta\Gamma^a{}_{bc}=\frac12g^{ad}
 (\nabla_bh_{cd}+\nabla_ch_{bd}-\nabla_dh_{bc}),
\qquad
\delta R^a{}_{bcd}=\nabla_c\delta\Gamma^a{}_{bd}
 -\nabla_d\delta\Gamma^a{}_{bc}.
\]

降低首指标还必须加入 \(h_{ae}R^e{}_{bcd}\)。在 Kinnersley 标架的
极端 Weyl 投影中，Ricci 扣除项和一阶标架修正不贡献，因此直接计算
\(\delta R_{lmlm}\) 与 \(\delta R_{n\bar m n\bar m}\)。
轨道 r0=6、a=.6、ell=m=2，在 r=4.5 和8、theta=1.1，恢复输入曲率
的相对误差小于6.4e-13。任意纯规范张量也经过零极端曲率检验。

## trace-driven kappa 的 resolvent 推导

引入辅助质量平方参数 lambda，保持原物理 Fourier 频率固定，解

\[
(\Box-\lambda)h_\lambda=16\pi T.
\]

对 lambda 求导并令 lambda=0，得到

\[
\Box\left.\partial_\lambda h_\lambda\right|_0=h_0,
\qquad
\kappa=-i\omega\left.\partial_\lambda h_\lambda\right|_0,
\qquad \Box\kappa=-i\omega h_0.
\]

使用对参数也求导的 retarded 边界条件；不能只对局部微分方程求导而
保持错误的无穷远波数。径向 Green 算子允许正负辅助质量平方，
角算符同时采用 a^2(omega^2-lambda)。这是辅助解析延拓，不是给物理云
引入负质量平方。单个 ell 的源投影依赖 lambda；额外接触项只有完整
角向求和才抵消。因此本验证只证明远离点源的方程，不能替代跳跃条件。

四阶差分为 \(D_4=(4D_{\delta/2}-D_\delta)/3\)，其中
\(D_\delta=[h_\delta-h_{-\delta}]/(2\delta)\)。检查 delta=1e-4、5e-5、
2.5e-5，并改变积分容差及外边界1000/2000，场值和梯度相对差小于4.3e-9。
这是参数稳健性检查，尚未显示单纯 delta^4 收敛率；差异已受数值积分误差影响。

## 复解与尚未确定的源归一化

AAB 构造首先给出复解。论文 Sec. III A 明确允许最后取实部。
对 Fourier 系数，实部操作是

\[
h^{\rm physical}_m=\tfrac12(h^{\rm complex}_m+
 \overline{h^{\rm complex}_{-m}}).
\]

不能直接把单个 Fourier 系数替换为其数值实部。也不能要求中间复解
满足真实场的正负 m 共轭关系。原始共轭差报告仅作为诊断保留。

独立张量计算确认了完整 Maxwell 场
F_ab=2 nabla_[a xi_b] 的圆周恒等式：

\[
\mathcal L_T\xi_a=2\nabla^b(f_{[a}{}^cF_{|c|b]})
+f_a{}^b j_b+\nabla_a(\xi_t-\tfrac12 f^{bc}F_{bc}),
\qquad j_a=\nabla^bF_{ba}.
\]

单凭这个恒等式，曾额外乘入自对偶系数2，使实现的总系数成为4。
现新增独立反算：对重构出的 xi 求 F=d xi，再收缩得到
F_lm 和 F_mbar n，分别与输入 phi0、phi2 比较。系数4在视界侧和
无穷远侧均给出精确的两倍，因此将总系数修正为2。两个投影恢复输入
的相对差小于3e-14，见 `spin1_circularity_development.json`。
这项检验能够约束齐次规范向量的幅度；真空 Einstein 和 Lorenz 检验不能。

下一项必需验收是完整角向展开下的点源连续性、径向导数跳跃及
相应 Einstein 分布源。还需静态 m=0、低多极和质量/角动量补全。
目前禁止将此开发合成结果视为已复现论文的环境源或通量。

初步匹配诊断 `metric_matching_development.json` 使用 theta=1.1、
r=6±5e-4，对物理实部 Fourier 系数求和到 ell=2..6。
两侧张量差相对于最大分量从约0.89降至0.54，但最大分量同时从约153
增长到2328，不能据此声称收敛。这一有限截断序列既未验证连续性，
也未验证正确的分布源；需进一步区分角向高阶接触项与辅助标量源组合。

## 角向投影进一步定位问题

`trace_radial_audit.json` 比较两个独立径向实现的 trace Green 核：
ell=2,4,6,8,10，r=4.5,6,8，相对差小于9.2e-12。Wronskian 的
变化也在同一量级。因此该范围内的径向积分不是已观察到的大差异的解释。

另对任意光滑规范向量独立构造 F=d xi、j_a=nabla^b F_ba，验证
O_(±1) phi_(0,2)=S_(0,2) j。该有源 Maxwell 检验通过，不能通过
任意调整源算子 S 的倍数来修复合成度规。

`metric_projected_matching_q8.json` 用8点 Gauss 积分把物理 Fourier
系数 h_tt 投影到 P_2^2(cos theta)。到 ell=5 时，两侧差分别为：
spin2 约+0.197456，spin1 约+0.879324，spin0 约-1.956105。
总差约-0.879325，几乎恰好少一个 spin1 分量。
这是定位线索，不是允许用匹配数据拟合系数的依据。

因此，上述把系数4改为2的反算检验仅证明**自对偶 Maxwell 分量**
归一化正确，不能推出完整实度规的规范向量已经正确。
必须继续推导反自对偶分量及复源的组合。当前合成度规仍为未验收开发实现。
