# 从局部点质量场独立推导标量模态的ell^-2尾部

这是对当前数值场的独立局部诊断，不是整篇复现验收。推导采用
Lorenz规范、测地线点粒子、光滑真空背景度规与光滑的有
质量背景标量Phi0。令m_p表示小天体质量，避免与方位模态m混淆。

## 1. Lorenz奇异度规的主导项

在粒子瞬时静止的局部惯性坐标中，主导线性Einstein方程为

\[
\nabla^2\bar h_{ab}=-16\pi m_p u_a u_b\delta^{(3)}(X),
\qquad \nabla^2(1/s)=-4\pi\delta^{(3)}(X).
\]

因此bar h_ab=4m_p u_a u_b/s。由于u^a u_a=-1，迹反转给出

\[
h_{ab}^{S}=\frac{2m_p}{s}(g_{ab}+2u_a u_b)+O(s^0).
\]

s是与世界线正交的局部空间距离。这是标准Lorenz点质量场的
主导局部形式；相关背景见Poisson、Pound、Vega的综述
[The Motion of Point Particles in Curved Spacetime](https://arxiv.org/abs/1102.0529)。
下面的标量尾部推导和数值检查是在当前项目中独立进行的。

## 2. 收缩到有源Klein–Gordon方程

记H_ab=∇a∇b Phi0，并使用背景方程g^ab H_ab=mu^2 Phi0。
文章的受迫源在粒子附近成为

\[
h_S^{ab}H_{ab}=\frac{2m_p}{s}A+O(s^0),\qquad
A=\mu^2\Phi_0+2u^a u^bH_{ab}.
\]

主导奇异阶上，局部空间Laplacian占优。因为∇²s=2/s，一个
对应的局部特解为

\[
\delta\Phi^S=m_p A s+\text{更高局部阶及光滑项}.
\]

质量项-mu²deltaPhi和沿世界线的有限频率时间导数都是较低
奇异阶。这个特解是连续但有尖点的场，而不是1/s发散的场。

## 3. 圆轨道上的A可以独立计算

u^a=u^t(1,0,0,Omega_p)，沿测地线有

\[
u^a u^b\nabla_a\nabla_b\Phi_0
=\frac{d^2\Phi_0}{d\tau_p^2}
=-(u^t)^2(\omega_c-\Omega_p)^2\Phi_0,
\]

其中背景云m_b=1。故

\[
A=[\mu^2-2(u^t)^2(\omega_c-\Omega_p)^2]\Phi_0.
\]

脚本同时通过协变Hessian直接收缩核对这个圆轨道表达式，并
检查u·u=-1以及局部奇异度规的收缩系数2A。

## 4. 转为球谐高阶尾部

在固定BL时间、r=r_p的赤道附近，以gamma表示单位球面上
离粒子的角距离，beta表示局部方向。主导距离为

\[
s=\gamma q(\beta)+O(\gamma^2),\quad
q(\beta)^2=g_{\theta\theta}\cos^2\beta+
(g_{\varphi\varphi}+u_\varphi^2)\sin^2\beta.
\]

球谐加法定理的核只依赖gamma，因此主导角向平均系数是
qbar=(2pi)^-1 integral q(beta) d beta。为了显式计算这个尖点的
尾部，可用具有相同主导尖点的模型2qbar sin(gamma/2)：

\[
\frac{2\ell+1}{2}\int_{-1}^1
\sqrt{2(1-x)}P_\ell(x)\,dx
=-\frac{4}{(2\ell-1)(2\ell+3)}.
\]

这个积分也可由Rodrigues公式和Beta积分得到：

\[
\int_{-1}^1(1-x)^\nu P_\ell(x)dx
=\frac{(-1)^\ell2^{\nu+1}\Gamma(\nu+1)^2}
{\Gamma(\nu-\ell+1)\Gamma(\nu+\ell+2)},
\]

将结果延拓至nu=1/2，再乘sqrt(2)(2ell+1)/2，即得到上式。
令L=ell+1/2，可得主导**球谐**模态预测

\[
\delta\Phi_\ell\sim-\frac{m_p A\bar q}{L^2}.
\]

这里不是把完整局部场等同于上述弦长模型，也没有据此固定
所有次领先项。脚本用独立积分核对Legendre系数，并用完全
椭圆积分核对qbar。

## 5. rp=20M的数值检查

单位Killing云质量、alpha=.3时，Phi0=-.00393154213815，
u^t=1.08286701398582，A=.000396241338996，qbar=20.2849441765。
按epsilon q归一后的独立预测系数为

\[
C_{\rm local}=A\bar q/\alpha^3=0.297693831184.
\]

对现有完整偶数阶场，定义C_ell=-L² deltaPhi_ell/(epsilon q)：

| ell | Re(C_ell) | Im(C_ell) | 与局部预测的复相对差 |
|---|---:|---:|---:|
| 6 | .361222918 | .010044789 | 21.61% |
| 8 | .329705297 | .000815945 | 10.76% |
| 10 | .316842717 | .000047519 | 6.43% |
| 12 | .310552822 | .000002079 | 4.32% |

数值序列趋向该独立预测的方向是有用的一致性证据，但不能视为
已经达到渐近极限。尤其当前数据是频率依赖的**椭球谐**分阶和，
推导使用球谐核；两种分阶方式的渐近对应尚未完成证明。有限
度规截断、源分辨率和次领先局部项也仍存在。

这没有消除与论文图7的差异，也没有引入任何拟合归一化。
可复现计算为report_local_scalar_cusp.py，数据见local_scalar_cusp.json。

后续已将完整偶数有限场（含ell'=0）直接重投影为球谐分阶，
结果见SCALAR_SPHERICAL_REPROJECTION.md。重建相对差约1e-14，
L=12的球谐缩放系数为.3105523307，与上表相近。该步骤核对
了有限输入的角基底转换；无限模态极限和度规截断仍未验证。
