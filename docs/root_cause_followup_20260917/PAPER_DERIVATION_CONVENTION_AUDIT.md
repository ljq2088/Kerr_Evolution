# 论文推导—实现审计：协变源、印刷遗漏、频率与规范

2026-09-17。本轮不修改生产物理代码、不拟合振幅或相位、不删除任何物理模态。所有本地路径相对仓库根目录。

**最明确的新发现是 Li v2 公开式 (66c) 与其定义式 (27) 的一个可解析定位的差异。现有生产实现直接计算正确的 `h^{ab}∇a∇bΦ`，没有这项遗漏。** 本轮通过 90 个独立度规分量输入验证差异，并量化了故意删除该项的反事实影响。不能据此断言作者真实代码执行了印刷表达式，也不能把该差异当成目前所有残差的原因。

证据类别：**原文明确**＝固定版本文本/公开源码；**推导等价**＝从几何定义推导并独立检验；**未知**＝公开材料不足以锁定实际作者运行流程。以下均区分这三类。

## 1. 固定来源与可复核位置

- Dyson 原论文 [arXiv:2501.09806v1](https://arxiv.org/abs/2501.09806v1)，源文件 `outputs/paper_original_reference/arxiv_v1_source/main_PRL.tex`（下称 D）。
- Li 后续独立工作 [arXiv:2507.02045v2](https://arxiv.org/html/2507.02045v2)，源文件 `outputs/paper_original_reference/li_2507_02045v2/source/main.tex`（下称 L）。源附录是 **Appendix B**，式 (65)、(66a–c)；Chandrasekhar 算子是 Appendix A。
- [Wardell–Kavanagh–Dolan, arXiv:2406.12510](https://arxiv.org/abs/2406.12510)，本地 `outputs/lorenz_reference/DirectLorenzGauge.tex`（下称 W）；其本地内容散列固定在 manifest，不仅依赖可变网页。
- [Dolan 等圆轨道重构 arXiv:2306.16459](https://arxiv.org/abs/2306.16459)，本地 `outputs/lorenz_reference/LorenzGaugeKerrCirc.tex`（下称 C）。
- [Sam 原公开生成器](https://github.com/srd24/KerrLorenzCirc/tree/5c1b42793ff893fd0c65a5b48ceecc02d34c5e35)，固定 commit `5c1b42793ff893fd0c65a5b48ceecc02d34c5e35`，`outputs/paper_metric_reference/srd24_KerrLorenzCirc/metric_reconstruction_calc_radiative.m`（下称 S）。这不证明 Li/Dyson 当时运行的文件逐字等于该快照。

## 2. δ□、源符号与迹反转

取 `δg_ab=h_ab`、`bar h_ab=h_ab−g_ab h/2`，固定标量 Φ。由连接变分直接得

\[
\delta g^{ab}=-h^{ab},\qquad
 g^{ab}\delta\Gamma^c_{ab}=\nabla_a\bar h^{ac},
\]
\[
\delta[(\Box-\mu^2)\Phi]
=-h^{ab}\nabla_a\nabla_b\Phi-(\nabla_a\bar h^{ac})\nabla_c\Phi.
\]

故线性扰动满足

\[
(\Box-\mu^2)\Phi_1=h^{ab}\nabla_a\nabla_b\Phi_0
 +(\nabla_a\bar h^{ac})\nabla_c\Phi_0.
\]

Lorenz 条件消去第二项。**收缩使用 h，不是 bar h。** 这是 L `eq:EOM_scalar_11_source` 179–188、Lorenz 349–352、式 (27) 371–375，与 D 式 (19)、`eq:EOM_sourced` 307–309 的相同内容。代码 `src/environment_source.py:147–153` 正是双次升指标后与协变 Hessian 收缩，符号正确；156–172 乘 Σ 再作角向投影。

误传迹反转张量的局部影响可明确写成

\[
\bar h^{ab}H_{ab}-h^{ab}H_{ab}=-\frac12h\Box\Phi_0
=-\frac12\mu^2h\Phi_0
\]

（末式仅在背景精确 on shell 时成立）。这是一个可单独投影的加性源，通常随 r、θ 和模态而变，并不产生统一的约 21.6 倍因子。当前 API 明确要求 covariant BL h，未在桥接接口中作迹反转。

**Dyson 的 off-shell 印刷不一致。** D:279 明确定义非密度算子 `Q=(□−μ²)Φ`；D:292–294 的式 (17) 除上述两项外又写 `−h Q/2`。对其显示的非密度定义直接变分没有该额外项；若重新定义密度算子，必须连同算子定义和其他变分一起改变，不能只额外补这一项。该项在精确云上为零，因此不改变 Kerr 的式 (19)。冻结虚频后 Q≠0 时，不能以此印刷式为依据随意加源。D:272 显示的作用量动能/质量项相对号也不能直接变分得到其 D:279 的 KG 符号；本实现以显示的物理 KG 方程、分离径向方程及边界定义为准，不为追逐这个印刷符号去改方程。

## 3. Li 式 (66c) 的逐分量恒等式检查

### 3.1 定义固定，避免升降指标和 Γ 权重混淆

L 式 (25)/(26)，355–369 给 normalized/unnormalized tetrad；式 (28) 与其后权重，378–396 给实际文件输出。令

\[
\Sigma=\Gamma\bar\Gamma,
\quad X\equiv\Sigma\Delta h_{l_+l_-},
\quad H\equiv g^{ab}h_{ab}.
\]

十个独立输入按顺序取

\[
(h_{l_+l_+},h_{l_-l_-},h_{m_+m_+},h_{m_-m_-},
 \Gamma h_{l_+m_+},\bar\Gamma h_{l_+m_-},
 \bar\Gamma h_{l_-m_+},\Gamma h_{l_-m_-},X,H).
\]

L:394–396 要求 `h_m+m−=ΣH−Δh_l+l−`。测试一次只把一个上述输入设为 1，其余为 0，利用显式 tetrad 矩阵逆变换回 covariant BL h；再分别计算 BL 双升指标收缩和正确使用 tetrad 度量升指标的收缩。两者最大绝对差 `1.18e−15`。

这是一项**点态张量恒等式测试**。输入不被声称为完整点粒子解。Lorenz 条件约束度规的一阶导数，并不限制任一点十个 h 分量的代数值；因此检验展开 `h:H` 的系数时，无需把每个单位张量误作全局物理解。

L 正文423把 L、L† 分别称为沿 m+、m−；但其式 (64)，829–837 中 `δ=L†/(sqrt(2)Γ)` 要求相反的角向记号。以显式 tetrad 为准，本测试的正确对照口径是

\[
\mathcal D_0=\partial_r-iK_c/\Delta,\quad
\mathcal D_0^\dagger=\partial_r+iK_c/\Delta,
\quad K_c=(r^2+a^2)\omega_c-am_c,
\]
\[
\mathcal L_0=\partial_\theta+Q_c,\quad
\mathcal L_0^\dagger=\partial_\theta-Q_c,
\quad Q_c=m_c\csc\theta-a\omega_c\sin\theta.
\]

这里算子作用于**云频率和云 m**，不能误用驱动度规或响应场的频率。程序同时保留正文423字面的另一种口径，未通过拟合挑选记号。式 (64) 的 NP Δ 首项还少了显式 tetrad (25) 所要求的 1/2；(63g,h) 的幂指数/共轭亦有显式记号问题。它们提示应锚定几何定义，而非从印刷表逐字修改生产。

### 3.2 遗漏项的解析推导

设 `F=R_c(r)S_c(θ)e^(−iω_ct+im_cφ)`，`B_ab=∇a∇bF` 是未归一化 tetrad 投影。其 tetrad 度量满足

\[
\eta_{+-}=2\Sigma/\Delta,\qquad\eta_{m_+m_-}=2\Sigma.
\]

仅考虑 X 和 H 两个输入，**正确的**源贡献为

\[
\Sigma\,h^{ab}B_{ab}\big|_{X,H}
=\frac{H}{2}B_{m_+m_-}
 +\frac{X}{2\Sigma^2}
   (\Delta B_{l_+l_-}-B_{m_+m_-}). \tag{R1}
\]

定义

\[
\mathcal E_r R=\partial_r(\Delta\partial_rR)+K_c^2R/\Delta,
\quad
\mathcal E_\theta S=S''+\cot\theta S'-Q_c^2S.
\]

直接从连接或 tetrad 求导可得

\[
 B_{m_+m_-}=\frac{2r\Delta}{\Sigma}R'_cS_c
 +\frac{2a^2\sin\theta\cos\theta}{\Sigma}R_cS'_c
 +R_c\mathcal E_\theta S_c,
\]
\[
 \Box F=(\mathcal E_rR_c\,S_c+R_c\mathcal E_\theta S_c)/\Sigma.
\]

将其代入 R1，X 的系数为

\[
\frac{\mathcal E_rR_c\,S_c-R_c\mathcal E_\theta S_c}{2\Sigma^2}
-\frac{2r\Delta R'_cS_c+2a^2\sin\theta\cos\theta R_cS'_c}{\Sigma^3}.\tag{R2}
\]

**印刷式 (66c)，L:889–912，恰好缺少 R2 的第一项 `X E_r[R_c]S_c/(2Σ²)`，其余项一致。** 全相位因子按式 (65) 乘回。对精确背景 KG，分离径向方程给

\[
\mathcal E_rR_c=(\mu^2r^2+\lambda_c)R_c,
\quad\lambda_c=A_c+a^2\omega_c^2-2am_c\omega_c.
\]

因此需要补在 `Σ S_Φ` 上的项是

\[
\boxed{\delta(\Sigma\mathcal S_\Phi)
 =\frac{X(\mu^2r^2+\lambda_c)}{2\Sigma^2}\Phi_c
 =\frac{\Delta h_{l_+l_-}(\mu^2r^2+\lambda_c)}{2\Sigma}\Phi_c.}\tag{R3}
\]

这**不正比** KG 残差，不能通过声明背景 on shell 消掉。若云故意 off shell，使用 R2 中的 `E_rR_c` 版本，不能不加检查地套用 R3。

### 3.3 无数值积分的解析反例及回归

取 a=μ=ω=m=0，Schwarzschild 的精确静态无质量解

\[
\Phi=(r-1)\cos\theta,\qquad \Box\Phi=0,
\]

令 X=1，trace 及其余八个独立分量为0。写 `f=1−2/r`、`H_0=1/(fr⁴)`，相应 BL 分量是

\[
h_{tt}=-f^2H_0/2,\quad h_{rr}=H_0/2,\quad
h_{\theta\theta}=-\Delta H_0/2,\quad
h_{\phi\phi}=\sin^2\theta h_{\theta\theta}.
\]

直接得到

\[
\Sigma h:H=2\cos\theta/r^4,
\qquad (\Sigma\mathcal S)_{\rm printed}=(3-r)\cos\theta/r^4.
\]

在 r=3、θ=.71，正确值 `0.01872498459235822`，印刷结果0。差就是 `(r−1)cosθ/r⁴`，与 R3 的 λ=2 一致。这也证明该印刷缺项不是 Kerr 独有或 massive 独有。

`audit_li_appendix_unit_basis.py` 完成4个任意光滑制造函数点、3个生产精确同步 Kerr 云点、2个上述解析 Schwarzschild 点，共90个单位分量检查。九个未涉 X 的通道误差、以及补回解析项后的全部通道误差，按每点十系数无穷范数归一，最大 `7.35e−16`。精确 Kerr 云的 R3 与一般 R2 遗漏项相差不超过 `1.70e−21`。脚本带 `<=1e−12` 恒等式断言、精确 KG 断言和上述闭式反例断言；它不以接近零系数为分母夸大误差。

## 4. 实际源和响应：只删该项的反事实实验

脚本 `audit_li_appendix_omission_response.py` 读取旧 Lmetric=6、角向12、径向 nr8/horizon32 的 **96节点**真实度规缓存。两轨道缓存均全部命中，无新度规求解。重新投影正确源与旧 J 的相对 L2 差分别 `2.61e−12`、`3.12e−12`。用 R3 投影得到 δJ，并比较

\[
J_{\rm correct},\quad J_{\rm printed\ counterfactual}=J_{\rm correct}-\delta J.
\]

共用原来的 Green、物理归一化、源范围和连续积分器；生产源不变。场幅度采用原 `R22 S22(pi/2)/α³`，不作任何校准。

| 轨道 r_p/M | Z_H 复振幅相对改变 | H通量变后/原 | I通量变后/原 | r=100场振幅变后/原 |
|---:|---:|---:|---:|---:|
| 42.1 | 30.66% | 0.5035 | 两者均0（束缚） | 1.6340 |
| 41.1 | 30.63% | 0.5041 | 2.7902 | 1.7011 |

r_p=42.1 的径向振幅谷由 `51.6803,141.0781` 变成 `15.6391,41.1551,141.1555`；远域的第二个谷几乎不移，但复场近似翻转符号，整体幅度增加。线性叠加检查相对误差≤`3.94e−10`，Green Wronskian变化≤`3.15e−11`。

这不是小修正或统一的21.6倍缩放，且不能把其影响说成作者数值结果的真实错误。详见 `li_appendix_omission_rp42p1.json`、`li_appendix_omission_rp41p1.json` 和对照图 `li_appendix_omission_control.png`。

### 4.1 与Li原图相同像素的单22替换

继续复用原18个场模态，**只把22换成上述反事实响应**；其余17个模态、各复相位、时间0、alpha^-3规范和Li颜色区间均不动。另按原来的正6选择重复。定义 `D_RMS=|| |Φ_local|−Φ_Li,nominal ||₂ / ||Φ_Li,nominal||₂`；全局数值使用预先固定9半径×72角，单环按72角。结果如下：

| 轨道 | 模态集合 | 原D_RMS | 只改22后D_RMS |
|---:|---|---:|---:|
| 41.1 | 全18 | 0.221383 | 0.719454 |
| 41.1 | 正6 | 0.126653 | 0.757659 |
| 42.1 | 全18 | 0.383290 | 0.794891 |
| 42.1 | 正6 | 0.346987 | 0.863190 |

全18集合的固定观测圆环细节：

| r_p | 观测r | 原D_RMS | 只改22后D_RMS |
|---:|---:|---:|---:|
| 41.1 | 50 | 0.157955 | 0.904964 |
| 41.1 | 100 | 0.057797 | 0.862894 |
| 41.1 | 150 | 0.139056 | 0.976610 |
| 42.1 | 50 | 0.631438 | 0.511434 |
| 42.1 | 100 | 0.155232 | 0.956032 |
| 42.1 | 150 | 0.595314 | 0.219232 |

个别圆环确实改善，但两个轨道/两个集合的全局误差均明显增大，不能挑改善的圆环作因果证据。r_p41.1总∞轨道能流由 `3.26663029e−5` 变为 `5.74771579e−5`（约+75.95%）；r_p42.1该总∞流保持 `1.83425591e−5`，因为22束缚，不向∞输送能量。

**在其余量固定的明确范围内，“仅此22印刷缺项解释当前场图偏差”的假说不受数据支持。** 这不等于已检验完整印刷公式对所有18个模态的影响，也不能证明作者具体用了哪个表达式。图 `li_appendix_single22_field_compare.png` 的灰带保留原色条量化区间，没有改振幅、旋转方位角或重新选择模式。

## 5. χ、κ、spin-0/1 与非静态偶极

| 公式/文本 → 本地实现 | 审计结论与边界 |
|---|---|
| W `eq:h-Lor` 369–370：M∂t hL=hAAB−2∇(ξ)；`lorenz_metric.py:74,113,133` | Fourier `e^(−iωt)` 给逆时间导数 `i/(Mω)`；M=1 的代码系数一致。静态不可套此除法。 |
| W 392、397：ξtrace=f·dh/2+dκ，□κ=M∂t h；`lorenz_metric.py:123–129`、`lorenz_kappa.py:3–6,43–55` | Resolvent `(□−λ)hλ=16πT` 求导给 □∂λhλ=h，所以 κ=−iω∂λhλ 正确。接触项与截断下源投影随λ变化另须检查，代码已有说明。 |
| C `eq:Box-kappa` 993–999：另一κ满足□κ=h/2 | 同名不等于同变量。离轨道真空区 `κ_W−2M∂tκ_C` 为齐次解，匹配条件可决定额外部分。不能仅以两文件中κ差2或ω就判错。 |
| W `eq:chi-DKW` 459–460；代码69 | ∂t²=−ω²，故代码的 `−1/(18ω²)` 与该定义一致。 |
| W 565–567、594–595：梯度中出现 χDKWSE+χSE−κ；代码128–133 | 代码用 `+∂(κ−χcompact)` 再作 `−2∇(ξ)`，相对号一致。χcompact源振幅还涉及分部积分/adjoint算子，不能由这一结构检查宣布所有跳跃系数已验证。 |
| W 537、593 与 `lorenz_metric.py:92–114` | antisymmetric two-form、divergence索引方向、完整两手性必须一起对照；不能只看代码`−2`就判倍数错误。当前模式启用`full_current=True`，两手性独立源；本轮没有重复其已完成的原作者张量数值核验。 |

**非静态 metric ell_g=1 的原始证据。** S:237–240 明确定义 s0下界|m|、s1下界max(1,|m|)、s2下界max(2,|m|)。540/542、868、890–891、1209、1375 分别在 spin1向量、匹配未知量、|m|=1闭合、径向求解和实际组装保留 ell1。L:377、402、785 明确使用这类原Notebook；因此没有“Lorenz重构仅ell_g≥2”的文本依据。D/L场图的ell≥2指**标量响应**，不可用于删metric偶极。作者环境代码的精确运行版本仍未知。

残余规范 Ward 恒等式（Ricci平直、精确云）为

\[
h^\xi_{ab}=-2\nabla_{(a}\xi_{b)},\quad
(\Box-\mu^2)(-\xi\cdot d\Phi_0)
=h_\xi:H-(\Box\xi)\cdot d\Phi_0.
\]

若□ξ=0且边界/轨道拼接均合法，局部源响应是纯规范场；但粒子内外分别纯规范不自动等于全域光滑齐次ξ。必须检验轨道处δ、δ′项及边界衰减，才能推出Z∞应为零。只删偶极比对图片不是规范证明，也不是旧公开门限bug的完整复刻。

## 6. Imω：原文确述、实际代码与未知顺序

- **Dyson Kerr明确**：D:262、689 假定同步阈值纯|211>，径向云就在阈值频率求解；693 用 KG Hessian迹检验。生产 `ThresholdCloud` 65–96 按此求真实同步a、实ω云。这一情形不存在需手动冻结的增长率。
- **Dyson Schwarzschild明确**：D:382（并在610重复）把Imω手设0，且承认不再是精确KG背景、渐近规范/守恒受影响；697独立在ingoing EF核验。**未知**：这句话未逐项规定是在径向谱解前取实，还是只冻时间。
- **本地Schwarzschild明确**：`environment_schwarzschild_cloud.py:68–95` 先保留复本征值求径向云和质量规范，仅`self.omega`取实以冻结时间；`report_environment_forced_mode.py:161–165` 记录这一选择。这是透明的实现约定，不能声称文本锁定了此唯一顺序。
- **Li明确**：L:312、320–326、332 用Leaver、a=.88、显示实数近似ωc≈.296294；L:729与Dyson a=.877的区别明确。全文没有明确说明其Imω在径向/时间/规范哪个阶段去掉。`environment_kerr_quasibound_cloud.py:8–9,75–87,124–134` 的新控制显式记录“先复云/质量、后仅冻时间”，属于对未知约定的控制，尚非作者流程证据。
- **旧阈值前凹陷不是复现目标**：D:374 已明确说波长在阈值附近极大、需远处提取，原图阈值前小凹陷为数值artifact。不能为贴这个凹陷而故意缩小远边界。该结论不能自动推广为Li全部场残差都是其边界误差。

## 7. 场、相干和、坐标与归一化

**R而非u。** L `eq:phi_11_decomp` 438–448 定义Φ11=ΣRSmode；541–555的 `hat R=sqrt(r²+a²)R` 只用于tortoise方程，其源为ΔJ/(r²+a²)^(3/2)，后面的Green式回到R。D:468/480亦以1/r渐近径向场分解。场图声称|Φ11|，没有u或|Φ|²定义。将sqrt(r²+a²)错乘进场/源会产生强r依赖，不会给3.5..200间近常数21.6。

**相干复和。** `environment_wake.py:69–77` 先加 `R S exp[i(mφ−ωt)]`，外部再取abs。实metric的负mg系数可共轭正mg，但云本身固定mc=1为复场；负mg响应频率是ωc−mgΩp，不是正响应的负频率，所以不能共轭整条标量通道。L:626正文说图中保留能贡献∞通量的模式，而图注明示更广的m和；当前同时保存18允许模态和6正模态，不能按哪组更贴图再决定。

频率选择ω=ωc+(m−1)Ωp给

\[
\Phi=e^{-i(\omega_c-\Omega_p)t}F(r,\theta,\phi-\Omega_pt).
\]

对实ω，改变t仅旋转每个圆环，不改圆环max/RMS。若一致变换 `t'=t+f(r), φ'=φ+g(r)`，径向系数变 `R'_m=e^(iω_m f−im g)R_m`；同样是每半径共同相位加方位角平移，不能解释圆环RMS的常数尺度差。混用BL云与另一坐标metric而不变相位才是错误；现有chart审计已对一致变换作了独立验证。

**21.6不是已知归一化修复系数。** 既有 `docs/LI_SOURCE22_33_CONVENTION_AUDIT_20260917.md` 的独立矢量源图读数显示 |J|/α³ 与公开中间源图差约21.5–21.9，但两模态径向形状和独立云profile相符。L:164–170 的Mc与Φ10记号、源图中间规范未结合原作图脚本完全闭合。raw Leaver规范的逆幅度约21.41958只是一项可比较的已定义值，既不精确等于差值，也没有证据证明作者图用该规范；不得拟合应用。

另一文本不统一是 D:454印刷stress只有L:132–136的一半，而D:475–476模态flux含2ωk、D:515 current无half，与简单C=1/2读法不一致。若严格一致地更换stress常数C，在固定Mc下云振幅按C^(−1/2)改、flux系数按C改，物理flux抵消；不能单边把场乘sqrt2，更不可能从中严格得到21.6。原云质量r²与Σ测度差已量化约1.49e−4（α=.3），也不够大。

## 8. 对主任务新增独立几何检查的只读复审

本节审查的是主任务生成的代码及已保存结果，不冒称本代理另跑了第二份度规数值实验。

`src/report_metric_vacuum_ricci_audit.py` 中 `ddG[q,a,b,c]` 表示 `∂q δΓa_bc`，33–37 行由 `δg^{-1}=−g^{-1}hg^{-1}` 和原始度规偏导直接构建。38–40 行六项恰为

\[
\delta R_{bc}=\partial_a\delta\Gamma^a_{bc}
-\partial_c\delta\Gamma^a_{ba}
+\delta\Gamma^a_{ad}\Gamma^d_{bc}
+\Gamma^a_{ad}\delta\Gamma^d_{bc}
-\delta\Gamma^a_{cd}\Gamma^d_{ba}
-\Gamma^a_{cd}\delta\Gamma^d_{ba}.
\]

二阶偏导的排列、连接乘积指标、逆度规导数及47行Lorenz散度升指标/迹项均逐项检查一致。`E δR E^T` 是tensor投影，应无复共轭；后续复Euclidean norm只作为固定normalized-null-tetrad基下的诊断尺度，不是Lorentz不变量。55–57行纯规范fixture为 `+Lξg`，因背景Ricci=0，其δR必须为0；不需要ξ满足□ξ=0，这与检查Lorenz纯规范的条件不同。

保存结果 `metric_vacuum_ricci_audit.json`：实际L6度规在r=3,20,50的δR/逐项范数和依次为 `7.98e−16,1.05e−13,3.60e−12`；r20将Jet order10降到8给 `1.17e−13`。纯规范fixture为`1.37e−16`。

`src/report_kg_metric_variation_audit.py` 从完整改变的度规 `g±εh` 和连接重算非密度KG算子，再作中心差分；其52–54行 `−h:H−C·dΦ` 的指标和符号正确，使用的是明确非Lorenz、off-shell制造输入。保存结果 `kg_metric_variation_audit.json` 四点最佳有限差分误差均≤`8.70e−11`。这直接支持第2节的无额外`−hQ/2`变分式，不依赖将KG残差置0。

这些检查强烈约束**局部真空算子、张量微分和KG变分**。任意尺度的齐次度规仍可满足真空方程，纯规范齐次项也可通过；因此它们不证明全局粒子δ源的归一化、轨道jump或作者选择的残余规范。两检查的脚本/输出哈希已纳入manifest，供与后续版本区分。

## 9. 证据边界与下一步判据

1. **已定位并证明**：Li显示式(66c)少R3项；当前直接协变实现正确。公开式的错误不是当前代码的新缺陷，单项反事实不能冒称作者真实程序。
2. **已排除简单解释**：一个统一21.6经验因子、R↔u、仅全局时间/方位角平移、把scalarell下界移到metricell、把不同论文κ同名当同变量。
3. **仍需具体对照**：作者源数组的复相位/原始中间归一化、Li fixed a=.88的Imω处理顺序、残余规范的实际轨道匹配与静态completion。当前公开证据不能将它们判定为某方错误。
4. 单个标量22反事实替换后的同像素比较另存 `li_appendix_single22_field_compare.json`，只改变该通道，保留其余17通道、原颜色区间与固定坐标；这项结果只检验“只这个遗漏解释现有场差”的假说。

所有新增脚本只生成诊断产物。精确原文/代码散列、运行结果散列及证据矩阵见本目录 `paper_derivation_audit_manifest.json`。
