"""Build the Chinese derivation PDF. Requires reportlab, svglib, matplotlib.
Run from repository root: .venv/bin/python docs/build_derivation_pdf.py
"""
import io
import os
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from svglib.svglib import svg2rlg
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak, Flowable
from reportlab.graphics import renderPDF
from reportlab.lib.pagesizes import A4

ROOT = Path(__file__).resolve().parents[1]
FONT = os.environ.get('KERR_CJK_FONT', '/mnt/c/Windows/Fonts/simsun.ttc')
pdfmetrics.registerFont(TTFont('Chinese', FONT, subfontIndex=0))
body = ParagraphStyle('Body', fontName='Chinese', fontSize=10.6, leading=18, spaceAfter=9, wordWrap='CJK')
heading = ParagraphStyle('Heading', parent=body, fontSize=18, leading=26, textColor=colors.HexColor('#153e64'), spaceAfter=17)
small = ParagraphStyle('Small', parent=body, fontSize=9, leading=15)
story = []
counter = 0


class Equation(Flowable):
    def __init__(self, source, number):
        super().__init__()
        source = source.replace(r'\Box_g', r'\nabla_\mu\nabla^\mu').replace(r'\phantom{0=}', r'\quad')
        source = re.sub(r'\\mathbf ([a-zA-Z])', lambda m: r'\mathbf{'+m[1]+'}', source)
        buff = io.BytesIO()
        math_to_image('$'+source+'$', buff, prop=FontProperties(size=12), format='svg')
        buff.seek(0)
        self.drawing = svg2rlg(buff)
        self.factor = min(1., 430/self.drawing.width)
        self.width = 475
        self.height = self.drawing.height*self.factor+17
        self.number = number

    def draw(self):
        self.canv.saveState()
        self.canv.translate(4,8)
        self.canv.scale(self.factor,self.factor)
        renderPDF.draw(self.drawing,self.canv,0,0)
        self.canv.restoreState()
        self.canv.setFont('Times-Roman',9)
        self.canv.drawRightString(480,self.height/2-2,f'({self.number})')


def p(text, style=body):
    story.append(Paragraph(text,style))


def eq(source):
    global counter
    counter += 1
    story.append(Equation(source,counter))


def page(title):
    if story:
        story.append(PageBreak())
    p(title,heading)


page('1　Kerr 标量场：从几何到波动方程')
p('双曲切片演化的逐步推导 · 2026-09-10',small)
p('本文对应项目 kerr-hyperboloidal 的实际实现，讨论固定 Kerr 背景上的无质量、无源、最小耦合标量测试场。采用 G=c=1、度规号差 (-,+,+,+)，限制 M&gt;0、|a|&lt;M。场的自旋权重为零；a 是黑洞旋转参数，与场的自旋权重不同。')
eq(r'\Box_g\Phi=\frac{1}{\sqrt{-g}}\partial_\mu(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\Phi)=0')
eq(r'\Delta=r^2-2Mr+a^2,\quad\Sigma=r^2+a^2\cos^2\theta,\quad L=r_+=M+\sqrt{M^2-a^2}')
p('从 Boyer-Lindquist 坐标引入 advanced ingoing 时间 v 与正则方位角 varphi。以下每次径向微分均须区分保持哪一个时间坐标不变。')
eq(r'dv=dt_{\rm BL}+\frac{r^2+a^2}{\Delta}dr,\qquad d\varphi=d\phi_{\rm BL}+\frac{a}{\Delta}dr')
p('在坐标顺序 (v,r,theta,varphi) 下，乘以 Sigma 的逆度规及体积因子为：')
eq(r'\Sigma g^{vv}=a^2\sin^2\theta,\quad\Sigma g^{vr}=r^2+a^2,\quad\Sigma g^{v\varphi}=a')
eq(r'\Sigma g^{rr}=\Delta,\quad\Sigma g^{r\varphi}=a,\quad\Sigma g^{\theta\theta}=1,\quad\Sigma g^{\varphi\varphi}=\csc^2\theta')
eq(r'\sqrt{-g}=\Sigma\sin\theta,\quad \Delta_{S^2}=\frac{1}{\sin\theta}\partial_\theta(\sin\theta\partial_\theta)+\frac{\partial_\varphi^2}{\sin^2\theta}')
p('将这些分量代入散度形式。交叉二阶导数出现两次；一阶项来自对 r²+a² 和 Delta 的径向求导：')
eq(r'0=a^2\sin^2\theta\,\Phi_{vv}+2(r^2+a^2)\Phi_{vr}+2a\Phi_{v\varphi}')
eq(r'\phantom{0=}+\Delta\Phi_{rr}+2a\Phi_{r\varphi}+2r\Phi_v+\Delta^{\prime}\Phi_r+\Delta_{S^2}\Phi')
p('这就是后续变换的起点。它含完整 Kerr 旋转效应，不是单独替换 Schwarzschild 有效势。式号是本文编号，与参考文献编号不同。',small)

page('2　场重标度：为什么演化 u=rΦ')
p('在未来零无穷附近，出射无质量场的主导行为通常为 Phi~F/r。定义 u=r Phi 后，辐射场 u 在紧致化端点可以保持有限。视界 r=L 不为零，因此这一变换也不引入视界奇点。')
eq(r'\Phi=\frac{u}{r},\quad\Phi_v=\frac{u_v}{r},\quad\Phi_r=\frac{u_r}{r}-\frac{u}{r^2}')
eq(r'\Phi_{rr}=\frac{u_{rr}}{r}-\frac{2u_r}{r^2}+\frac{2u}{r^3},\qquad\Phi_{vr}=\frac{u_{vr}}{r}-\frac{u_v}{r^2}')
eq(r'\Phi_{r\varphi}=\frac{u_{r\varphi}}{r}-\frac{u_\varphi}{r^2},\quad\Phi_{vv}=\frac{u_{vv}}{r},\quad\Delta_{S^2}\Phi=\frac{\Delta_{S^2}u}{r}')
p('代入上一页波动方程，再整体乘 r。径向一阶项、时间一阶项、零阶项分别由以下组合得到：')
eq(r'B=\Delta^{\prime}-\frac{2\Delta}{r}=2M-\frac{2a^2}{r}')
eq(r'2r-\frac{2(r^2+a^2)}{r}=-\frac{2a^2}{r},\qquad W=\frac{2\Delta}{r^2}-\frac{\Delta^{\prime}}{r}=-\frac{2M}{r}+\frac{2a^2}{r^2}')
p('因此正则场的方程为：')
eq(r'0=\Delta u_{rr}+2(r^2+a^2)u_{vr}+a^2\sin^2\theta\,u_{vv}+B u_r')
eq(r'\phantom{0=}+2a u_{r\varphi}+2a u_{v\varphi}-\frac{2a^2}{r}u_v-\frac{2a}{r}u_\varphi+W u+\Delta_{S^2}u')
p('这里的 u_r 仍是保持 v 不变的导数。下一步变换到双曲时间时，不能把它直接当作新切片上的径向导数。尤其 -2a/r 的方位角一阶项不能遗漏；固定 m 后它成为一个纯虚系数。')
p('本文保留有量纲几何时间 tau。u 是重标度后的标量场，不是引力波应变 h、Weyl 标量 Psi4 或能量通量。',small)

page('3　双曲时间与紧致化的链式法则')
p('取 horizon-penetrating hyperboloidal minimal gauge，将未来零无穷固定在 sigma=0，将未来事件视界固定在 sigma=1：')
eq(r'\tau=v-h(r),\quad h(r)=2r+4M\ln(r/L),\quad\sigma=\frac{L}{r}')
eq(r'H=h^{\prime}=2+\frac{4M}{r},\quad H^{\prime}=-\frac{4M}{r^2},\quad q=\frac{d\sigma}{dr}=-\frac{\sigma^2}{L}')
p('H 是高度函数 h 的径向导数，不是 h 本身。对任意 u(tau,sigma)，有：')
eq(r'\left.\partial_v\right|_r=\partial_\tau,\qquad\left.\partial_r\right|_v=q\partial_\sigma-H\partial_\tau\equiv J')
eq(r'u_{vr}=q u_{\tau\sigma}-H u_{\tau\tau},\qquad u_{r\varphi}=q u_{\sigma\varphi}-H u_{\tau\varphi}')
p('二阶径向导数需要对 q、H 再求导。注意 Jq=q q_sigma，JH=H′：')
eq(r'J^2u=q^2u_{\sigma\sigma}-2qH u_{\tau\sigma}+H^2u_{\tau\tau}+q q_{,\sigma}u_\sigma-H^{\prime}u_\tau')
eq(r'q_{,\sigma}=-\frac{2\sigma}{L},\qquad q q_{,\sigma}=\frac{2\sigma^3}{L^2}')
p('其中最后两项是变系数链式法则的贡献。若只将 J 当成常系数算子平方，将得到错误的一阶导数系数。')
p('与参考文献的记号对应：Ripley 的 ingoing 时间 v_R=v-r，故他的高度函数导数是 -1-4M/r，而这里相对于 advanced v 使用 -2-4M/r；二者给出同一个 tau，至多相差常数。他使用 rho=1/r，本项目使用 sigma=L rho。')
p('Macedo 的 radial function fixing minimal gauge 采用无量纲时间，关系为 tau_Macedo=tau/L。方位角使用 ingoing varphi，在整个变换中保持不变。',small)

page('4　整理系数：消去端点的表观奇异性')
p('将上一页导数代入 u 的 ingoing 方程，先按微分阶数收集，不立即展开 Delta：')
eq(r'A_0=\Delta H^2-2(r^2+a^2)H,\quad C=2q(r^2+a^2-\Delta H),\quad D=\Delta q^2')
eq(r'E_0=-\Delta H^{\prime}-BH-\frac{2a^2}{r},\qquad F_0=\Delta q q_{,\sigma}+Bq')
eq(r'0=(A_0+a^2\sin^2\theta)u_{\tau\tau}+C u_{\tau\sigma}+D u_{\sigma\sigma}+E_0u_\tau+F_0u_\sigma')
eq(r'\phantom{0=}+2a(1-H)u_{\tau\varphi}+2a q u_{\sigma\varphi}-\frac{2a}{r}u_\varphi+Wu+\Delta_{S^2}u')
p('直接在 sigma=0 计算 r=L/sigma 会产生无穷大相消。应先用下列恒等式化简，再在端点计算多项式：')
eq(r'\Delta H-2(r^2+a^2)=-8M^2+\frac{4Ma^2}{r}')
eq(r'A_0=H(-8M^2+4Ma^2\sigma/L),\qquad H=2+4M\sigma/L')
eq(r'C=2L+\frac{(2a^2-16M^2)\sigma^2}{L}+\frac{8Ma^2\sigma^3}{L^2}')
eq(r'D=\sigma^2\left(1-\frac{2M\sigma}{L}+\frac{a^2\sigma^2}{L^2}\right)')
eq(r'E_0=\frac{(2a^2-16M^2)\sigma}{L}+\frac{12Ma^2\sigma^2}{L^2}')
eq(r'F_0=2\sigma-\frac{6M\sigma^2}{L}+\frac{4a^2\sigma^3}{L^2},\qquad W=-\frac{2M\sigma}{L}+\frac{2a^2\sigma^2}{L^2}')
p('A0 的量纲为长度平方，C、E0 的量纲为长度，D、F0、W 无量纲，因而各项量纲一致。这些多项式就是代码中的径向系数。',small)

page('5　球谐投影与 Kerr 模态耦合')
p('固定整数 m，采用归一化球谐及 exp(im varphi) 约定。角向投影使用单位球面的坐标测度 dOmega=sin(theta)dtheta dvarphi，而不是 Kerr 截面的物理面积测度。')
eq(r'u=\sum_{\ell=|m|}^{\ell_{\max}}u_{\ell m}(\tau,\sigma)Y_{\ell m},\quad\Delta_{S^2}Y_{\ell m}=-\ell(\ell+1)Y_{\ell m},\quad\partial_\varphi Y_{\ell m}=imY_{\ell m}')
eq(r'S_{\ell\ell^{\prime}}=\int Y_{\ell m}^*\sin^2\theta\,Y_{\ell^{\prime}m}\,d\Omega,\qquad A=A_0 I+a^2 S')
eq(r'E=E_0+2iam(1-H),\quad F=F_0-\frac{2iam\sigma^2}{L}')
eq(r'V_{\ell\ell^{\prime}}=\left[W-\frac{2iam\sigma}{L}-\ell(\ell+1)\right]\delta_{\ell\ell^{\prime}}')
eq(r'A\ddot{\mathbf u}+C\partial_\sigma\dot{\mathbf u}+D\partial_\sigma^2\mathbf u+E\dot{\mathbf u}+F\partial_\sigma\mathbf u+V\mathbf u=0')
p('角向矩阵可以不经过数值积分而精确构造。令 X 为乘以 cos(theta) 的球谐矩阵：')
eq(r'\cos\theta\,Y_{\ell m}=c_{\ell+1}Y_{\ell+1,m}+c_\ell Y_{\ell-1,m},\quad c_\ell=\sqrt{\frac{\ell^2-m^2}{4\ell^2-1}}')
eq(r'S_{\ell\ell}=1-c_\ell^2-c_{\ell+1}^2,\quad S_{\ell,\ell+2}=-c_{\ell+1}c_{\ell+2},\quad S_{\ell,\ell-2}=-c_\ell c_{\ell-1}')
p('在最低模态处令 c_|m|=0。因此只有同一奇偶性的 l 相互耦合，m 不变。虽然 A 只含 l±2 的带状项，求 A 的逆后，加速度可依赖同一奇偶扇区的多个模态。')
p('计算 S=P(I-X²)P 时，应先包含中间模态 lmax+1，再取目标子矩阵。若先截断 X 再平方，最高模态的对角元会缺失 c_(lmax+1)²。代码对此用独立 Gauss-Legendre 积分检查。')
p('对照 src/kerr_scalar.py：ainv 保存逐径向点的 A 逆矩阵；C、D、E、F 是径向列向量；V 还含球谐本征值。整体乘 -1 并令 sigma=L rho，可逐项恢复 Ripley 式 (11) 的自旋权重零情形。',small)

page('6　类空切片、特征边界与零自旋极限')
p('先验证演化时间的几何性质。在连续角向上，时间二阶导数的系数为 A(theta)=A0+a²sin²(theta)。设 kappa=a/L，则 2M/L=1+kappa²，有：')
eq(r'\frac{A(\theta)}{L^2}=-4(1+\kappa^2)[1+(1+\kappa^2)\sigma][1+\kappa^2(1-\sigma)]+\kappa^2\sin^2\theta<0')
p('这证明 tau=常数切片在计算域类空。球谐 Galerkin 投影同样保持负定性，代码初始化时也检查 A 的全部本征值。')
eq(r'D=\sigma^2(1-\sigma)(1-\kappa^2\sigma),\quad D(0)=D(1)=0')
p('对径向主部代入局部平面波或 u=f(sigma-c tau)，特征速度满足：')
eq(r'A c^2-Cc+D=0,\qquad D=0:\quad c_1=0,\quad c_2=C/A')
eq(r'C(0)=2L>0,\qquad C(1)=-\frac{2(L^2+a^2)}{L}<0')
p('因此 sigma=0 处非零特征速度为负，流出 sigma≥0 的计算域；sigma=1 处非零特征速度为正，流出 sigma≤1 的计算域。另一个特征在边界相切。没有特征从端点进入外部演化域，因而不额外给定 Dirichlet、Neumann 或 Sommerfeld 数据。')
p('数值上仍在两个端点计算同一个演化方程，只是径向二阶项系数自动为零。端点场值并不固定为零；波形输出正是这些动态端点值。')
p('令 a=0、L=2M，所有 m 相关项和角向耦合消失。可得到便于独立核对的 Schwarzschild 系数：')
eq(r'A=-16M^2(1+\sigma)I,\quad C=4M(1-2\sigma^2),\quad D=\sigma^2(1-\sigma)')
eq(r'E=-8M\sigma,\quad F=2\sigma-3\sigma^2,\quad V_{\ell\ell}=-\sigma-\ell(\ell+1)')
p('本项目的网格截止于事件视界，没有演化黑洞内部。这里的特征分析是对该初值演化问题的边界说明；仅要求频域函数有界，并不足以在任意函数空间中唯一识别 QNM。',small)

page('7　空间谱离散与时间 Runge-Kutta 演化')
p('径向采用 N+1 个 Chebyshev-Lobatto 配点，包含两个物理端点：')
eq(r'\sigma_j=\frac{1-\cos(j\pi/N)}{2},\quad j=0,\ldots,N')
eq(r'w_j=(-1)^j\eta_j,\quad\eta_0=\eta_N=\frac{1}{2},\quad\eta_j=1\ (0<j<N)')
eq(r'(D_1)_{ij}=\frac{w_j}{w_i(\sigma_i-\sigma_j)}\ (i\ne j),\quad(D_1)_{ii}=-\sum_{j\ne i}(D_1)_{ij},\quad D_2=D_1D_1')
p('令 U 的行对应径向点、列对应 l；P=U_tau。径向微分通过左乘 D1、D2 完成，角向 A 的逆在每个径向点分别作用。对每个径向行 i：')
eq(r'\dot{\mathbf U}_i=\mathbf P_i')
eq(r'\dot{\mathbf P}_i=-A_i^{-1}\left[C_i(D_1P)_i+D_i(D_2U)_i+E_i\mathbf P_i+F_i(D_1U)_i+V_i\mathbf U_i\right]')
p('将 U、P 展平并拼接为 y，就得到常系数线性半离散系统 y_dot=L_disc y。这是 method of lines：空间先离散，再对时间积分。这里 D_i 是 PDE 的标量系数，不是微分矩阵 D1。')
eq(r'k_s=f\!\left(t_n+c_s\delta t,\ y_n+\delta t\sum_{j<s}a_{sj}k_j\right),\quad y_{n+1}=y_n+\delta t\sum_s b_s k_s')
p('代码使用 scipy.integrate.solve_ivp(method="DOP853")，即八阶显式 Runge-Kutta 方法。时间步长由误差估计自适应调整，不是固定步长 RK4。默认 rtol=10^-9、atol=10^-11，另设 max_step=0.5M。')
p('默认 401 个输出时间由 t_eval 指定，它们是存储时刻，并非内部积分步长。DOP853 的稠密输出用于生成这些时刻的解。')
p('Chebyshev 点在端点附近聚集，最小间距为 O(N^-2)。显式积分仍受离散算子的稳定性限制，自适应误差控制不能替代空间与时间收敛检查。当前代码没有加入人工耗散或吸收层。')
p('复数数组用于固定 m 的方位角扇区。要构造一个实标量解，可取该复解的实部；若显式同时拼接 ±m 扇区，则必须保持相应球谐共轭约定与归一化。',small)

page('8　初始数据、提取量与验证依据')
p('当前初值仅激发指定 l0，角向其余模态为零；径向是端点包络修饰的 Gaussian，而非纯 Gaussian：')
eq(r'u_{\ell m}(0,\sigma)=\delta_{\ell\ell_0}\left[\frac{\sigma(1-\sigma)}{c(1-c)}\right]^2\exp\left[-\left(\frac{\sigma-c}{w}\right)^2\right],\quad p_{\ell m}(0,\sigma)=0')
p('默认 c=0.45、w=0.12、l0=2。p=0 只表示初始双曲时间导数为零；不应据此称为几何时间反演对称数据，也不是纯出射或纯 QNM 初态。初始端点场值为零不等于以后施加零边界条件。')
eq(r'u^{\mathscr{I}^+}_{\ell m}(\tau)=u_{\ell m}(\tau,0),\qquad \Phi^{\mathcal{H}^+}_{\ell m}(\tau)=u_{\ell m}(\tau,1)/L')
p('原波形图左右两面板都画 u 的系数，右图若需物理 Phi 则除以 L；未指定单独角向观察点。三维图重建 Re[sum_l u_lm Y_lm]/r，在固定双曲时间而非固定 BL 时间上作图。')
p('已执行的验证：8 项单元测试通过，涉及系数与 Ripley 原式的独立对应、角向求积、径向多项式微分、端点方向、零自旋解耦、m 共轭与耦合、质量缩放和参数检查。')
p('M=1、a/M=0.7、tau≤80M 的端点波形比较：N=48 与 64 的最大差为 2.46×10^-9（m=0）和 9.87×10^-10（m=2）；lmax=6 与 8 的共同模态最大差为 1.55×10^-8 和 7.88×10^-9。它们是离散结果间的差，不是严格误差上界。详见 docs/validation.md。')
p('文献与可追溯来源',heading)
p('[1] R. P. Macedo, A. Zenginoglu. Hyperboloidal Approach to Quasinormal Modes. Frontiers in Physics 12, 1497601 (2025). DOI: 10.3389/fphy.2024.1497601；arXiv:2409.11478v2。框架入口：第 3.2 节。',small)
p('[2] J. L. Ripley. Computing the quasinormal modes and eigenfunctions for the Teukolsky equation using horizon penetrating, hyperboloidally compactified coordinates. CQG 39, 145009 (2022). DOI: 10.1088/1361-6382/ac776d；arXiv:2202.03837v2。直接系数依据：第 2 节式 (5)-(11)。',small)
p('[3] R. P. Macedo. Hyperboloidal framework for the Kerr spacetime. CQG 37, 065019 (2020). DOI: 10.1088/1361-6382/ab6e3e；arXiv:1910.13452v2。规范对应：第 4.1 节与附录 C。',small)
p('本推导以 Klein-Gordon 散度式独立展开，并对照上述文献。对应代码：src/kerr_scalar.py；测试：tests/test_kerr_scalar.py。尚未加入有质量场、源、非线性反作用，亦未完成 QNM 拟合、极长时尾波或能量通量守恒验证。',small)


def footer(canvas,doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor('#d5e0ea'))
    canvas.line(52,48,A4[0]-52,48)
    canvas.setFont('Chinese',8)
    canvas.setFillColor(colors.HexColor('#536779'))
    canvas.drawString(52,34,'Kerr 双曲切片标量场 · 推导与实现')
    canvas.drawRightString(A4[0]-52,34,str(doc.page))
    canvas.restoreState()


pdf = ROOT/'docs/kerr_hyperboloidal_derivation_zh.pdf'
doc = SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=52,leftMargin=52,topMargin=48,bottomMargin=62,
                        title='Kerr 双曲切片标量场：逐步推导与数值实现',author='Kerr Evolution project')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
print(pdf)
