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

## 显式反自对偶源补全

现已实现 `lorenz_spin1_chiral.py`。对反自对偶通道采用相反频率及 m
的 Teukolsky 齐次解，再共轭其伴随核，但**不共轭原复电流的源算子**。
这样积分得到的是同一个复电流的另一种手性，而不是另一个共轭电流。
两种通道均保留完整源 j_SE-j_DKWSE，其中

\[
\langle V,j_{\rm DKWSE}\rangle
 =\frac{8\pi}{u^t}\frac{4i}{9\omega}
 u^a u^b(\mathcal A^\dagger V)_{ab}.
\]

这里 A 是论文显式 DKW 电流算符，沿用已实现的伴随形式。
自对偶通道的此项在真空齐次核积分中为零，复现旧的接触项捷径；
反自对偶通道中此项非零，不能省略。最终将两个二形式相加再取散度。
没有根据轨道匹配结果拟合或手动加倍现有分量。

初步8点角积分、ellmax=4 给出的 h_tt 投影两侧相对差为2.20e-4，
相较只保留一种手性的约0.72已显著减小。径向导数跳跃为-0.3804684。
由点粒子源直接推导预期值：

\[
[\partial_r h_{tt,m}] =-\frac{8}{u^t\Delta_0}
 (u_t^2+\tfrac12g_{tt})\frac{\delta(\theta-\pi/2)}{\sin\theta},
\]

与 P_2^2 投影后得 -0.3806134392。径向采样仍相隔0.001，且模态与
积分都有限截断；这一接近结果还需精度收敛和其他张量分量验证。

加入两种手性后，`metric_development.json` 也记录了未强制对称化的
完整复解正负 m 共轭差：r=4.5 为3.5e-12，r=8 为6.2e-13。
这一自然满足的关系是额外交叉证据，不取代分布源验证。
旧 `metric_reality_development.json` 和 `metric_projected_matching_q8.json`
保留的是仅含自对偶自旋1部分时的历史诊断，不能当作当前实现的输出。

12点角积分的 ellmax=4 结果与8点积分一致。增加到 ellmax=6，投影
两侧相对差为1.80e-5，导数跳跃为-0.380527486，与理论值相差8.60e-5。
两侧相距0.001，因此场值本来就存在随间距缩小的一阶变化，不能要求
这一有限距离差随 ellmax 单调降到零。必须分别收敛径向间距和角模态截断。
固定 ellmax=4，把间距缩小十倍后，导数跳跃误差从1.45e-4降到6.78e-5，
而场值差没有下降，表明此时角截断仍是限制。尚需在更高 ellmax 下
重复径向间距序列，不能把单个投影的这些结果称为完整源验收。

## 全分量 Taylor 精度与静态 trace

全张量径向外推曾发现 order=6 不足：DKW 标量的四阶算符加上两个
梯度消耗六阶，某些分量的剩余 Taylor 系数已不可靠。现在全张量匹配
采用 order=8，保留二阶径向导数；独立两侧直接求值验证了外推。
原 h_tt 检查的两个指标均为时间指标，消耗的径向阶数较少，不能将其
阶数需求推广到所有分量。失败的 order=6 全分量结果不作为物理证据。

开始单独实现零频部分：`lorenz_static_trace.py` 给出 m=omega=0 的
精确 trace Green 函数。令 b=sqrt(1-a^2)、x=(r-1)/b，方程成为

\[
\frac{d}{dx}\left[(x^2-1)\frac{dR}{dx}\right]-\ell(\ell+1)R=0.
\]

视界正则解为 P_l(x)，无穷远衰减解为 Q_l(x)，径向 Wronskian 为 -b。
因此点源 J_l=-16 pi Y_l0(pi/2)/u^t 的 retarded 静态解为

\[
R_l(r)=-\frac{J_l}{b}P_l(x_<)Q_l(x_>).
\]

已验证 a=0,.6,.9、ell=0,2,4 的真空方程、连续性和精确源跳跃，
并检验奇宇称源为零及单极远端1/r行为。这仅是静态 trace，尚不是
完整静态度规或质量/角动量补全，不能用非静态的除频率公式生成后者。

`metric_tensor_matching_q12.json` 覆盖10个独立分量和 P22/P32 两种
投影。使用可靠 Taylor 阶数后的最大单侧极限值差，在 ellmax=2,4,6
分别为1.260、0.202、0.0256；最大导数跳跃误差分别为5.105、0.747、0.0931。
这些是 M=1 的坐标分量绝对误差，不是统一物理范数。
非 tt 分量仍有较大的角截断误差，需要更高模态及角积分分辨率，
当前尚不通过全张量点源验收。
