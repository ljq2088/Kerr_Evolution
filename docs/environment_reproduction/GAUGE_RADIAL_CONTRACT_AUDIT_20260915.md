# s=0、±1 径向调用与约定专项审核

参数：M=1，a=0.8771530275949366，rp=20，Omega_p=0.011071760582072043。

本报告检查频率、自旋、分离常数、基解归一化和高阶Taylor递推的契约。没有改生产代码，也不把契约通过解释成总视界通量32.66%偏差已解决。

## 生产调用链及上一轮控制的边界

| 路径 | 基解/源计算 | 是否随pybhpt后端控制改变 |
|---|---|---|
| lorenz_metric._homogeneous_radial_data → homogeneous_field_jet | 当前场点的In或Up，s=0、±1、±2 | 是 |
| lorenz_spin1.spin1_amplitudes | 轨道处同一s的In、Up及Wronskian，构造自对偶源幅度 | 是 |
| lorenz_spin1_chiral.chiral_amplitudes | 包括(-m,-omega)和另一手性的源幅度 | 是 |
| lorenz_chi.chi_amplitudes | 轨道处s=0的两支径向及Wronskian，构造紧支撑辅助标量源 | 是 |
| lorenz_kappa.scalar_resolvent → environment_radial.RadialGreen | 独立DOP853的有质量/无质量标量In、Up | 否 |
| lorenz_kappa.trace_field_jet | msq=0的trace以及质量延拓后的trace | 否 |
| lorenz_kappa.kappa_jet | msq=±step、±step/2的四次差分，step=min(5e-5,.01omega²) | 否 |

生产trace/kappa默认使用rtol=1e-12、rmax=2000、视界偏移1e-4。`environment_trace_variation.TraceMassVariation` 是独立诊断实现，不由生产度规自动导入；只有显式调用其安装函数才替换kappa。此前88节点AUTO/TEUK源比较明确保持这些路径固定，不能用该比较排除trace/kappa共同问题。

`lorenz_trace.trace_amplitudes` 还提供另一条RadialGreen标量路径（默认rmax=4000、offset=1e-6），但`nonstatic_metric`不调用它；不要把这一诊断入口的参数当作生产trace参数。

## 径向方程和lambda约定

采用exp(-i omega t+i m phi)，定义Delta=r²-2r+a²、K=(r²+a²)omega-am。把Teukolsky的径向散度项展开，得到

Delta R''+(s+1)Delta'R'+[(K²-2is(r-1)K)/Delta+4is omega r-lambda_s]R=0。

生产`separated_jet`正是此式；其angular方程使用

A_s=lambda_s-a²omega²+2amomega。

因此pybhpt径向lambda已包含a²omega²-2amomega，不能再补一次该项。pybhpt的SWSH接口也采用相应lambda约定。项目的`environment_cloud.angular_eigenvalue`返回的是标量A，`lorenz_kappa.trace_field_jet`先转换为lambda再传给`separated_jet`，这一转换正确。4组实际参数的独立标量A→lambda残差最大2.84e-13；spin1的lambda_+1=lambda_-1-2也通过。

质量延拓时额外径向势为-msq r²；角向项为a²(omega²-msq)cos²theta。它仍使用a²omega²-2amomega把A转换为径向lambda，不能把这一转换中的omega²也替换成omega²-msq。

## 负频率共轭与spin1 TS归一化

由于K在(m,omega)同时变号时反号，径向方程给出同一s下

R_s(-m,-omega)=conj[R_s(m,omega)]，

这里保持同一In/Up的单位透射归一化。测量了ell=1、6、18的m=1以及ell=m=6，三种自旋、8个径向点：物理分段R和R'的最大相对残差5.98e-14，未发现负频率或手性径向共轭错误。

另外独立展开D=d/dr-iK/Delta，利用负自旋径向方程直接计算

D²R_-1=R_-1''-2i(K/Delta)R_-1'-(i(K/Delta)'+(K/Delta)²)R_-1。

对于单位边界透射幅度，预测D²R_-1=C_bc R_+1，其中

C_In=-4K_+²-2iK_+(r_+-r_-)，

C_Up=-B²/(4omega²)，B²=lambda_-1²+4amomega-4a²omega²。

这些常数由视界Frobenius幂和无穷远波支推导，未用数值拟合。实际物理分段最大残差1.59e-13，主导ell=m=1两支分别1.07e-14和4.89e-14。它检验了自旋变换的约定；由于库本身也使用TS，不能独立排除TS和径向积分共同的所有误差。

## 源幅度与场的基解归一化必须一起变化

对于R_In→c_In R_In、R_Up→c_Up R_Up，Wronskian乘c_In c_Up。线性源算子的外区幅度乘1/c_Up，内区幅度乘1/c_In，所以强迫物理场不变。这说明一致的整体边界相位或整体幅度通常会在度规源幅度和场的乘积中抵消；真正风险是两处使用不一致的归一化，或基解混入错误波支。

本轮给每个s、正负m、In/Up施加彼此不同的任意复常数，同时覆盖源幅度和齐次场，并清除相关内存缓存。在r=10和40、theta=1.1的实际ell=m=1度规上：

| 扇区 | r=10相对度规范数变化 | r=40相对变化 |
|---|---:|---:|
| spin1，含两个手性和DKW源 | 1.89e-12 | 2.02e-14 |
| spin0，trace+kappa+chi | 6.29e-14 | 1.37e-14 |

未发现源幅度与齐次场的归一化、复共轭或缓存混用错误。spin0这项只改变pybhpt chi基解，trace/kappa原样保留，不能据此验证trace/kappa的物理源公式。

## 高阶径向Taylor递推

另用80位运算实现多项式形式

Delta²R''+(s+1)Delta Delta'R'+[K²-2is(r-1)K+Delta(4is omega r-lambda)]R=0

的卷积递推，未调用Jet，也未对生产的有理函数形式做逐阶微分。固定相同的双精度R、R'输入，检查第2至8阶导数。生产double Taylor的最大相对差1.30e-8，发生于s=+1、ell=6、m=1、r=19.9的第8阶导数；extended系数运算后整个样本最大5.32e-12。

对当前生产order=6实际能提供的第2至6阶，ell=1样本最大1.29e-10，ell=6为4.94e-12，ell=18为1.69e-15。小导数的相对误差会放大，且负自旋检查包含撤去Kinnersley的zeta因子，不能把这些数直接当作最终度规或通量误差。此前发现的度规大项相消还会进一步放大输入误差，需在重构层单独处理。

## 仍需独立排查的共同误差

1. trace/kappa的标量resolvent、质量导数、质量依赖角向投影及延拓边界条件；本轮pybhpt契约对这些没有覆盖。
2. HBL与TEUK共同采用的边界波支或同一渐近展开错误；TS/共轭契约可同时通过，需MST锚定的独立波支/Green核检验。
3. 紧支撑源的Lorenz伴随算子、两个手性的物理装配和delta源匹配。基解重标度不变不能证明这些共同公式正确。

复现：`PYTHONPATH=src OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_gauge_radial_contract_audit.py`。

原始数据：`gauge_radial_contract_audit.json`，包含源文件哈希、逐点残差和预测TS常数。
