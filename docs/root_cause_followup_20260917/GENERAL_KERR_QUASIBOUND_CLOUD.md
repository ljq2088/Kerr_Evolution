# 固定 Kerr 自旋的复频准束缚云：推导、接口与局部验证

日期：2026-09-17。范围：新增独立模块 `src/environment_kerr_quasibound_cloud.py`；不替换旧同步云，不修改既有源、Green 核或通量模块。本文验证云与局部源接口，未进行新的完整度规、源网格或环境场计算。

模块 SHA256：`3403dbf846cf630f6668fa2ed698272a37c687124a3d193d68b155b43816458b`。原始数值证据见同目录 `general_kerr_quasibound_cloud_validation.json`，由 `src/report_kerr_quasibound_cloud_validation.py` 生成；该 JSON 保留运行时依赖哈希，不随其他模块随后变动而改写。

## 1. 物理目标与频率

采用 $M=1$、$\mu=\alpha=0.3$、$a=0.88$、$\ell=m=1$，求解

$$
\Phi=e^{-i\omega_s t+im\phi}R(r)S(\theta),\qquad (\Box-\mu^2)\Phi=0.
$$

默认频率为此前独立 Leaver 连分式与双端复射击核验的复准束缚根：

$$
\omega_s=0.296293534725611462814885674586
+2.2166093964457949473\times10^{-9}i.
$$

本模块消费该谱值并检验双端匹配，不重新求谱。根的证据见 `docs/environment_reproduction/fixed_spin_cloud_frequency_20260917.json` 和 `fixed_spin_kerr_cloud_shooting_20260917.json`。默认根仅适用于上述参数；其他参数必须明确传入已独立确立的复频。目前仅对 $|211\rangle$ 提供经过验证的接口。

这与旧的 $a=0.8771530275949366$、$\omega=m\Omega_H$ 同步实频云是不同背景。不能仅修改受迫 Green 核中的 $a,\omega$ 就宣称已经使用固定 $a=0.88$ 的云。

## 2. 角向方程与两种内积

令 $c^2=a^2(\omega_s^2-\mu^2)$，使用分离常数 $A$ 的约定：

$$
{1\over\sin\theta}\partial_\theta(\sin\theta\,\partial_\theta S)
+\left[c^2\cos^2\theta-{m^2\over\sin^2\theta}+A\right]S=0.
$$

以标准归一化球谐在 $\phi=0$ 的函数 $Y_{\ell m}(\theta,0)$ 展开，矩阵为

$$
H_{\ell j}=\ell(\ell+1)\delta_{\ell j}-c^2(\cos^2\theta)_{\ell j}.
$$

复 $c^2$ 时这是复对称矩阵，不能当作 Hermitian 矩阵调用 `eigh`。新模块使用一般复本征问题，并固定目标球谐系数为正实数，保证相位可复现。

**能量归一化**需要 Hermitian 内积：

$$
2\pi\int_{-1}^{1}|S(z)|^2dz=\sum_j|b_j|^2=1.
$$

**复频角向模态投影**则来自复对称 Sturm–Liouville 算子的双线性正交关系：

$$
2\pi\int_{-1}^{1}S_j(z)S_k(z)dz=0\quad(j\ne k),\qquad
\widetilde S_j={S_j\over\sum_n b_{jn}^2}.
$$

方位角 Fourier 投影已经移除 $e^{im\phi}$ 后，极角源投影使用 $\widetilde S_j$，不是 $S_j^*$。两者在实频实角函数极限相同。`project_complex_source` 显式实现此双线性投影；它并未使既有实频 `RadialGreen` 或 `mode_flux` 自动支持复频。

验证用人为放大的 $c^2=0.1+0.2i$ 展示区别：正确双线性异模投影约 $10^{-15}$，错误 Hermitian 异模投影约 $-0.00860i$。这不是实际云的误差量；实际 $c^2=-0.001711533407079335+1.017200780987525\times10^{-9}i$。

## 3. 从 BL 径向方程到正则穿越视界方程

记

$$
\Delta=(r-r_+)(r-r_-),\quad r_\pm=1\pm\sqrt{1-a^2},\quad
K=(r^2+a^2)\omega_s-am.
$$

BL 径向方程为

$$
(\Delta R')'+\left[{K^2\over\Delta}-\mu^2r^2-a^2\omega_s^2+2am\omega_s-A\right]R=0.
$$

使用 ingoing Kerr–Schild 坐标

$$
T=t+h(r),\quad\Psi=\phi+g(r),\quad h'={2r\over\Delta},\quad g'={a\over\Delta},
$$

并固定积分常数

$$
h={2r_+\over d}\log{r-r_+\over2}-{2r_-\over d}\log{r-r_-\over2},\quad
g={a\over d}\log{r-r_+\over r-r_-},\quad d=r_+-r_-.
$$

此时 $r_*=r+h$，且

$$
F=e^{i\omega_s h-img}R,\qquad
\Phi=e^{-i\omega_s T+im\Psi}F(r)S(\theta).
$$

令 $H=(2r\omega_s-am)/\Delta$，代入

$$
R=e^{-i\omega_sh+img}F,\quad
R'=e^{-i\omega_sh+img}(F'-iHF)
$$

并收集 $F'',F',F$ 得到

$$
\boxed{\Delta F''+[\Delta'-2i(2r\omega_s-am)]F'
+[(\omega_s^2-\mu^2)r^2+2r\omega_s^2-A-2i\omega_s]F=0.}
$$

代码直接积分这个正则方程，避免先积分奇异 BL 导数、再通过大数相消恢复 $F'$。

未来视界入射条件 $R\sim e^{-i(\omega_s-m\Omega_H)r_*}$ 在上述坐标和相位常数下等价于

$$
F_H=\exp\left[-i\omega_s r_++im\left({a\over2}+{a\over r_+}\log{d\over2}\right)\right].
$$

这由坐标变换直接确定，没有拟合任意相位。$a=0$ 时化为 $e^{-2i\omega_s}$，与已有 Schwarzschild 云的相位一致。

令 $x=r-r_+$、$F=\sum c_nx^n$，定义

$$
b_0=d-2i(2r_+\omega_s-am),\quad b_1=2-4i\omega_s,
$$

$$
v_0=(\omega_s^2-\mu^2)r_+^2+2r_+\omega_s^2-A-2i\omega_s,
\quad v_1=2r_+(\omega_s^2-\mu^2)+2\omega_s^2,\quad
v_2=\omega_s^2-\mu^2.
$$

以 $c_0=F_H$、负下标系数为零，递推为

$$
c_n=-{[(n-1)(n-2)+b_1(n-1)+v_0]c_{n-1}+v_1c_{n-2}+v_2c_{n-3}
\over n[(n-1)d+b_0]}.
$$

无穷远采用衰减支 $k=i\sqrt{\mu^2-\omega_s^2}$、$\Re\sqrt{\mu^2-\omega_s^2}>0$，使用现有 massive scalar 逆半径渐近级数取得对数导数，再转换为 $F,F'$。两端用 DOP853 积到中间半径，以振幅连接并独立检验对数导数；匹配误差超过 $10^{-8}$ 时拒绝该谱值/边界组合。

## 4. 穿越未来视界切片上的质量归一化

采用复标量作用量对应的 $C=1$ 约定：

$$
T_{\mu\nu}=\partial_\mu\Phi^*\partial_\nu\Phi+
\partial_\mu\Phi\partial_\nu\Phi^*
-g_{\mu\nu}(|\nabla\Phi|^2+\mu^2|\Phi|^2),
$$

$$
j^\mu=-i(\Phi^*\nabla^\mu\Phi-\Phi\nabla^\mu\Phi^*).
$$

Kerr–Schild 坐标中 $\sqrt{-g}=\Sigma\sin\theta$、$\Sigma=r^2+a^2\cos^2\theta$，所需逆度规为

$$
g^{TT}=-{\Sigma+2r\over\Sigma},\quad g^{Tr}={2r\over\Sigma},\quad g^{T\Psi}=0,
\quad g^{rr}={\Delta\over\Sigma},\quad g^{r\Psi}={a\over\Sigma},
\quad g^{\theta\theta}={1\over\Sigma},\quad g^{\Psi\Psi}={1\over\Sigma\sin^2\theta}.
$$

在 $T=0$ 上定义外部云的 Killing 能量 $E=-\int\Sigma T^T{}_T\,dr\,dz\,d\Psi$ 和 Noether 电荷 $Q=\int\Sigma j^T\,dr\,dz\,d\Psi$，其中 $z=\cos\theta$。代入 $\partial_T\Phi=-i\omega_s\Phi$、$\partial_\Psi\Phi=im\Phi$，并记 $X=\Im(F^*F')$，得到

$$
E=2\pi\int dr\,dz\left\{
\left[(\Sigma+2r)|\omega_s|^2+\mu^2\Sigma+{m^2\over1-z^2}\right]|F|^2|S|^2
+\Delta|F'|^2|S|^2+|F|^2|S_\theta|^2+2amX|S|^2\right\},
$$

$$
Q=4\pi\int dr\,dz\left\{(\Sigma+2r)\Re\omega_s|F|^2+2rX\right\}|S|^2.
$$

先计算原始能量 $E_{\rm raw}$，再乘振幅 $\sqrt{M_c/E_{\rm raw}}$。没有使用 $E=\omega Q$ 的捷径，也没有沿 BL $t=\mathrm{const}$ 奇异切片归一化。数值积分区间为 $[r_++\epsilon,R_{\rm out}]$；有限截断和求积误差通过独立加密及通量平衡检验，不将它们宣称为严格零。

穿越视界的谱模态通量为

$$
F_E^H=2[2r_+|\omega_s|^2-am\Re\omega_s]|F_H^{\rm norm}|^2,
\quad F_Q^H=2[2r_+\Re\omega_s-am]|F_H^{\rm norm}|^2.
$$

指数增长/衰减因子要求

$$
2\Im\omega_s E+F_E^H=0,\qquad 2\Im\omega_s Q+F_Q^H=0.
$$

这给归一化、相位、复频和视界边界一个独立的整体校验。新的 KS 积分式还与直接构造完整逆度规和应力张量的局部计算逐点核对。

## 5. 明确的增长冻结近似

默认 `freeze_growth=False`：时间、径向与角向均使用复谱频率，是数值误差范围内的准束缚本征模。

`freeze_growth=True` 的顺序固定为：

1. 用复 $\omega_s$ 求 $R,S$ 及空间导数；
2. 用复谱模态在 KS $T=0$ 上的 Killing 能量归一化；
3. 仅将时间导数和受迫 Fourier 频率中的云频率设为 $\omega_R=\Re\omega_s$，保留复空间轮廓和上述振幅。

因此 `.spectral_omega` 总是复数，冻结时 `.omega` 为实数。这是一项显式近似，不是精确实频本征模，也不是同步云；在保持 BL 复空间轮廓的同时冻结时间，并不定义一个处处正则的精确未来视界本征解，只能在指定有限外部源区间按近似评估。

保持 BL 空间函数不变时，其 KG 缺陷可独立推导为

$$
{(\Box_{\omega_R}-\mu^2)\Phi\over\Phi}
=g^{tt}(\omega_s^2-\omega_R^2)+2g^{t\phi}m(\omega_R-\omega_s).
$$

新 Hessian 和 Jet 保留谱空间导数并仅改变时间导数，确实再现这个缺陷。缺陷会随近视界 BL 系数放大；不能只根据 $\Im\omega_s$ 很小就称残差为零。

`ks_radial()`、`ks_integrals()` 默认返回/归一化底层复谱模态；显式传入 `frequency=cloud.omega` 可检查实际冻结近似在同一切片上的量。代码不会再次隐藏地归一化冻结后的能量。作者原文允许何种冻结顺序需由文献约定另行核对；此模块不会将本选项自动标记为作者唯一算法。

## 6. 接口与现有源计算的连接

```python
from environment_kerr_quasibound_cloud import GeneralKerrQuasiboundCloud
cloud = GeneralKerrQuasiboundCloud(a=.88, alpha=.3, mass=1., freeze_growth=True)
R, Rp, Rpp = cloud.radial_state([3., 20., 100.])
field, hessian, inverse_metric = cloud.hessian(20., 1.1)
```

| 接口 | 含义 |
|---|---|
| `a,mu,ell,m,rp,rm,rmin,rmax,match` | 背景、模式与实际已解积分域 |
| `spectral_omega,omega,c2,lam` | 复谱频率、使用的时间频率、谱角参数与 $A$ |
| `mass,charge,amplitude,normalization` | 请求质量、谱电荷、振幅、归一化约定 |
| `radial`, `radial_state` | BL $R,R'$ 或 $R,R',R''$；复数组 |
| `ks_radial` | 正则 $F,F'$；可显式改变时间频率作诊断 |
| `angular`, `hessian`, `lorenz_source` | 复角函数、协变 Hessian、$h^{\mu\nu}\nabla_\mu\nabla_\nu\Phi$ |
| `jet(...,return_geometry=True)` | 谱空间 Taylor 系数与正确时间频率的 GHP 几何 |
| `ks_integrals`, `horizon_balance` | 独立能量、电荷和增长/视界通量平衡 |
| `provenance` | 模块哈希、谱/时间频率、切片、相位、冻结顺序与求解参数 |

默认 $r_{\min}=1.4749746834815165$、$r_{\max}=70/\Re\kappa_c\simeq1489$；所有 BL/KS 径向入口拒绝非有限、积分域外半径，避免 dense output 外推。极角接口限于离轴 $0<\theta<\pi$ 的源求积点。

冻结频率时，已有 `environment_source.project_source` 可通过新云的自定义 Hessian 回调使用本模块；两个半径上的固定对称 toy 度规接口测试与显式投影逐复数完全一致。它只验证调用约定，不是假装已经验证真实度规源。复频精确模式则须使用复角向双线性投影和将来相应复频受迫径向方法，不能直接套用既有实频通量。

任何新源文件必须纳入本模块的 `provenance`。旧通用 `source_fingerprint()` 并不自动包含本新增模块；不能将旧源重新标记为新云。父任务将使用单独的 `background_provenance`。

## 7. 已完成的数值验证

精化设置：$\epsilon=10^{-6}$，无穷远 $75/\Re\kappa_c$，DOP853 相对容差 $2\times10^{-13}$，视界六阶、无穷远八阶，归一化径向点 6001，角求积 48 点。粗设置同时取 $\epsilon=10^{-5}$、$60/\Re\kappa_c$、容差 $2\times10^{-11}$、四/六阶、3001 点；该联合变化用于稳定性核对，不单独归因于某一个参数。

| 检验 | 结果与适用范围 |
|---|---|
| 双端对数导数匹配 | $-2.99\times10^{-15}-1.73\times10^{-17}i$ |
| 6001→12001 点独立能量求积 | $1.0000000000000002\to1.000000000029166$ |
| 能量/电荷视界平衡相对残差 | $1.44\times10^{-13}$ / $2.09\times10^{-12}$ |
| 粗/精 $R,R',R''$ | 五个半径的最大相对差 $1.58\times10^{-10}$ |
| 独立 BL Frobenius 视界态与相位 | 相对差 $5.73\times10^{-15}$，无相位拟合 |
| 直接 KS 应力张量与推导积分密度 | 三点相对差 $\le4.62\times10^{-16}$ |
| 精确模式协变 KG 残差 / $\mu^2|\Phi|$ | $r_++0.0005$ 为 $9.95\times10^{-13}$；其余三点 $\le1.31\times10^{-15}$ |
| 五点差分独立径向 ODE 残差 | 以各 ODE 项绝对值和归一，最大 $1.03\times10^{-10}$ |
| Jet 与 Hessian 的 KG 运算 | 最大差 $1.38\times10^{-12}$，近视界点主导 |
| 同步极限对旧 `ThresholdCloud` | 无相位拟合，径向最大 $2.39\times10^{-10}$，Hessian 最大 $2.81\times10^{-10}$ |
| $a=0$ 对旧 `SchwarzschildCloud` | BL 径向最大 $1.05\times10^{-10}$；KS 导数近视界比较 $8.12\times10^{-8}$ 受旧 BL→KS 相消影响 |
| $a=0$ 新云能量通量平衡 | 相对残差 $1.28\times10^{-10}$ |
| 源接口 | 固定 toy 对称度规、两半径，新显式投影与旧接口冻结分支逐复数相等 |
| 有针对性的测试 | `pytest -q tests/test_kerr_quasibound_cloud.py`：10 passed in 0.73s |

冻结近似的 KG 缺陷 / $\mu^2|\Phi|$ 为：

| 半径 | 冻结缺陷 |
|---|---:|
| $r_++0.0005$ | $7.07104\times10^{-7}$ |
| 3 | $3.18304\times10^{-8}$ |
| 20 | $1.62035\times10^{-8}$ |
| 100 | $1.48926\times10^{-8}$ |

这些数值与上述解析缺陷相符，差值的归一化最大 $1.76\times10^{-12}$。冻结近似在同一 KS 切片上的能量为 $1.0000000269886269$，底层谱模态为 $1.000000000029166$；约 $2.70\times10^{-8}$ 的差别被记录，没有再吸收进振幅。

以上结果建立了新的云和局部导数接口的可用性。它们不等于完整环境场、通量、有限 $\ell$ 源投影或论文图像已经一致；实际 $a=0.88$ 度规与新云的受迫源仍需另做有 provenance 的计算。
