# Li 后续工作：统一对象、符号、归一化与方法

参考固定为 arXiv:2507.02045v2。状态：对齐实现进行中，尚未达到完整方法一致，暂不发布新的最终通量比较。不能将旧同步背景结果或 KS 归一化结果改标签视为 Li 同参数结果。

## 1. 对象与参数

- 几何单位 G=c=1，数值 M=1；号差 (-,+,+,+)；BL 坐标 (t,r,theta,phi)。
- 主黑洞固定 a/M=0.88，mu M=alpha=0.3；不再用同步根 a=0.877153... 代替。
- 背景分别为 (ell_c,m_c,n_c)=(1,1,0) 与 (2,2,0)，即通常的 |211> 与 |322>。论文 n_c 是径向泛音，不能与氢原子主量子数混淆。
- 顺行圆赤道轨道；Omega_g=1/(sqrt(r0^3/M)+a)。研究区间 r0/M=[4.1,50.1]，文章未给全体实际采样点。
- 偶极场图用 r0=41.1,42.1，t=0，赤道 theta=pi/2，r<=200M；图注 theta=0 与极角定义矛盾，按明确的赤道几何解释并保留记录。

## 2. 符号与微扰阶数

Li epsilon=mp/M，等于旧代码 q；Li zeta=alpha^3 sqrt(eta)，eta=Mc/M，等于旧 Dyson epsilon。不得只按变量名搬用旧 epsilon。

Phi=zeta Phi^(1,0)+epsilon*zeta Phi^(1,1)+...；g=gKerr+epsilon*h^(0,1)+zeta^2*h^(2,0)+...。

omega_g=m_g*Omega_g 是度规 Fourier 频率，Omega_g 是轨道频率。源选择规则 m=m_c+m_g，omega=omega_c+m_g*Omega_g。论文通量处 omega_g/Omega_g 的命名应以物理频移 omega-omega_c 为准，禁止再重复乘一个 m_g。

全局 Fourier 约定 exp(-i omega t+i m phi)。真实度规 H_-mg=conj(H_mg)，这里不添加2；背景云是复标量，负 m 的标量响应不能由正 m 响应直接共轭得到。云整体相位、t=0及phi=0须写入元数据，不能用逐模相位去贴图。

## 3. 角向与分离常数

2*pi*integral_0^pi |S_lm|^2 sin(theta) dtheta=1；球谐也使用相同单位球面范数及同一相位约定。

角向方程：(1/sin theta)d_theta(sin theta*d_theta S)+[a^2(omega^2-mu^2)cos^2 theta-m^2/sin^2 theta+Lambda]S=0。

径向常数 lambda=Lambda+a^2*omega^2-2*a*m*omega。当前旧函数参数 lam 指 Lambda；新记录必须用 angular_eigenvalue_Lambda 与 radial_separation_lambda 区分。

若保留复 omega，角向矩阵是复对称而非 Hermitian，不能直接套用实频角函数的 Hermitian 正交关系。这与作者实际忽略虚部的步骤有关；其精确处理顺序尚未公开，不能擅自称作一致。

## 4. 四标架与度规数据

Gamma=r+i a cos(theta)，barGamma=r-i a cos(theta)，Sigma=Gamma*barGamma，Delta=r^2-2Mr+a^2。

l_±=[±(r^2+a^2)/Delta,1,0,±a/Delta]，m_±=[±i a sin(theta),0,1,±i csc(theta)]。

正规化基：l=l_+，n=-Delta/(2Sigma)l_-，m=m_+/(sqrt(2)Gamma)，bar m=m_-/(sqrt(2)barGamma)，故 l.n=-1，m.bar m=1。

作者输出的加权分量顺序固定为：
1. h_l+l+；2. h_l-l-；3. h_m+m+；4. h_m-m-；
5. Gamma*h_l+m+；6. barGamma*h_l+m-；
7. barGamma*h_l-m+；8. Gamma*h_l-m-；
9. Sigma*Delta*h_l+l-；10. trace h。

h_m+m-=Sigma*h-Delta*h_l+l-。源使用升指标后的 h 而非迹反转 bar h。作者第10项若仍是扁球谐系数，必须先变换为球谐后，才能与前9项一并作球谐截断。

方法目标是 KerrLorenzCirc 输出加权自旋球谐 j_g<=18；独立重构的扁球谐 L<=18不自动等价。静态 mg=0及质量、自旋补全单独验收。

## 5. 云谱与质量归一化

作者用 Leaver 法，径向连分式/级数基准150项；加密参考另外记录。频率采用实际谱根，不把印刷约数0.296294当作完整精度输入。

应按论文 Eq(zeta_def)核对：Mc=-integral(T^t_t * Sigma dr dOmega)。应力能为
T_ab=partial_a Phi partial_b conj(Phi)+partial_b Phi partial_a conj(Phi)-g_ab(|grad Phi|^2+mu^2|Phi|^2)，没有实标量常见的额外1/2。

现有 GeneralKerrQuasiboundCloud 是复谱空间轮廓、KS T=0质量归一化，再可选冻结时间增长；它只作为独立对照，不能冒充上述 BL 归一化。复准束缚态的视界极限与先冻结/后冻结并非自动可交换，必须显式保存积分下限及其变化、频率处理和归一化证据。

## 6. 源的方法对齐

核心几何式 S=h^{ab}nabla_a nabla_b Phi_c；径向源是投影 Sigma*S，而不是 S本身。

Li 路线是先把加权四标架分量和云导数写成径向×角向项，再展开 Gamma^{-beta}barGamma^{-sigma}，预计算角耦合系数和径向 Fourier 系数。

Sigma^{-beta}=sum_{p=0}^{12}f_p(r)cos(p theta)，奇数p严格为0；f_p使用64位十进制运算。这不是整个求解器或浮点输入有64位物理精度。

Gamma分组后正幂展开必须保留binomial(n,k)。印刷二项式少写该系数，不作为删除它的理由。

附录混合迹分量的漏项已由独立坐标张量收缩证实。正确源、按印刷式的诊断源分别保存；不为接近参考图而删项，也不从公式差异推断作者实际程序必定有同样错误。

新 li_separable_factors.py 实现上述 Fourier 因子，是完整半解析源的基础模块。旧二维 Gauss 源仅保留为独立验证，不能在配置中改名为半解析源。

## 7. 径向与边界

R是物理标量径向函数；hat R=sqrt(r^2+a^2)R。dr*/dr=(r^2+a^2)/Delta；两种变量及其导数、Wronskian必须明确转换。

In=(r-r+)^[-i(2Mr+omega-a*m)/(r+-r-)] sum A_j(r-r+)^j，A0=1，j<=4。印刷式若缺m，应按径向方程和视界生成元得到omega-m*OmegaH。

Up=exp(i k r*) r^(i M mu^2/k)/sqrt(r^2+a^2) sum B_j/r^j，B0=1，j<=4。正频传播 k>0，负频传播按出射群速度选择符号；束缚频率选衰减支。精确k=0必须单独处理。

W0=Delta*(Rin*Rup'-Rup*Rin')恒定；R=(Rup*integral_inner^r Rin*J dr+Rin*integral_r^outer Rup*J dr)/W0。

端点位置、容差未完整公开：明确列作自主数值选择，并通过端点移动验收。不能把当前Coulomb/Whittaker方法直接称为作者四阶方法，也不能选不收敛但更贴图的边界。

## 8. 观测量、符号、单位

定义N为离开云的电荷流：无穷远向外为正、视界向内为正。N与论文表示云电荷变化的Qdot符号相反。

Ninf=2*sgn(omega)*sqrt(omega^2-mu^2)*|Zinf|^2（传播支）；NH=4Mr+*(omega-m*OmegaH)*|ZH|^2。

Fwave=omega*N；有效轨道能量损失Fscalar=(omega-omega_c)*N=m_g*Omega_g*N；轨道方程中 dEorb/dt=-sum(Fscalar+FGW)。所以负的有效视界通量意味着轨道从这一渠道获能，不等于普通波能量流也必为负。

Lscalar=(m-m_c)N=Fscalar/Omega_g。两者作为归一化/频率选择的独立恒等式检查。

单位云质量场hatPhi对应物理扰动epsilon*sqrt(eta)*hatPhi；论文系数Phi^(1,1)=hatPhi/alpha^3。单位epsilon^2*eta的通量除以alpha^6才是单位epsilon^2*zeta^2的论文曲线。复场先求和再取模；不同度规阶驱动同一标量模式时先加复振幅再平方。

总通量两端l<=5；rp20逐模图至l=10；场图l=2..5。场包含不向无穷远传播的束缚模式；不能按传播判据筛掉场模态。图注“全部模态”和正文正m表述存在歧义，主清单按全部允许模态，另列文字筛选诊断，不能择优贴图。

## 9. 执行与验收顺序

1. 完成本文件的可核实配置与未公开项登记。
2. Leaver背景、BL质量条件及复频处理对照。
3. 同参数作者球谐度规、完整半解析源、四阶边界的单通道闭合。
4. 补齐每个目标的模式，分别检查源、角向、径向、边界收敛。
5. 最后对比源图、逐模通量、总通量、场图；缺模不报总量。

li_alignment_contract.py只检验声明配置，不是物理正确性的自动证明。target_contract.json保存原文散列及待办。未解决项不得以空值通过门槛；此前结果不重新贴标。

## 本轮实际完成（2026-09-25）

- 新增 li_alignment_contract.py：声明配置检查、符号映射、频率选择和场/通量单位转换。旧自旋、二维源、KS质量归一化不能通过Li配置门槛。配置通过仍不代表实际算法正确。
- 新增 li_separable_factors.py：pmax=12、64位径向Fourier系数；以独立高精度角积分验算。
- 新增 li_separated_source.py 和 li_source_terms.json：86个印刷源展开项、3个由协变恒等式要求补回的单项式；预计算角耦合，以球谐度规径向系数组装J。完整十个分量使用独立坐标Hessian单元检验。
- 在原作者a=0.877153、rp20、j<=3保护模式中，新方法与直接源最大相对差1.8531e-13（pmax12）。这些是方法验证输入，绝不算作a=.88生产结果。
- 新增 li_leaver_cloud.py：对a=.88、alpha=.3的两种云态实际求Leaver谱及径向级数，保留复杂频率，不暗中冻结或质量重标度。150项偶极根Re(omega)=0.2962935347281123，四极根Re(omega)=0.2984513875640803；300项结果另存。所检半径2..320的有限级数ODE残差最大分别7.81e-10和4.12e-11（150项），不是全部空间严格误差界。

质量归一化、实际全角阶作者度规输入、四阶边界统一驱动和全部模式计算仍未完成，因此本轮没有发布新的最终通量对比。


后续执行更新：已新增有限 BL 归一化、统一四阶 Green 模块并启动完整通量批次。当前状态、具体数值与未核实项以 [FLUX_PROGRESS_20260925.md](FLUX_PROGRESS_20260925.md) 和批次 execution.json 为准；此前完成情况描述仅为当时快照。
