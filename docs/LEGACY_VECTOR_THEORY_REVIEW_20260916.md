# 旧版 vector 缺项诊断：理论复核（2026-09-16）

复核对象：用户提供的 `Downloads/Kerr_environment_horizon_flux_diagnosis.md`。本次只读核验原论文 TeX、固定版本公开代码、现有归一化报告，并读取其他审计模块完成的消融结果；本复核没有重跑 Mathematica 或通量积分，没有更改生产模块。

**结论：公开代码中漏掉非静态 \(\ell=1,|m_g|=1\) vector 扇区的门限错误已经有直接证据；但“它单独解释 2025 年论文视界偏差”已被本轮消融及下述符号界反证（§5）。** 报告开头、§12、§23 的确定语气需要收紧。特别需要分清真空区域的 Lorenz 约束、轨道处的匹配条件，以及复现旧计算链与给出合法物理修正这三件事。

## 1. 删 vector：真空区域仍可满足 Lorenz，但轨道处不能忽略匹配和分布项

令 \(\bar h^{ab}=h^{ab}-g^{ab}h/2\)，\(C^b[h]=\nabla_a\bar h^{ab}\)。在背景满足 \((\Box-\mu^2)\Phi_0=0\) 时，完整的一阶标量源为

\[
 (\Box-\mu^2)\delta\Phi
 =h^{ab}\nabla_a\nabla_b\Phi_0+C^b[h]\nabla_b\Phi_0.
\]

这是 2025 论文 `main_PRL.tex:290–309` 中约束化简前后的直接关系。生产函数 `environment_source.py:147–153` 只计算第一项，因此要求输入度规满足相应的 Lorenz 条件。

原重构论文明确把 spin-one 扇区定义为 **traceless vacuum Lorenz pure gauge**（`LorenzGaugeKerrCirc.tex:962–979`）：

\[
 h^V_{ab}=-2\nabla_{(a}\xi_{b)},\quad
 \nabla_a\xi^a=0,\quad \Box\xi_b=0
 \quad (R_{ab}=0,\ r\ne r_p).
\]

于是 \(h^V=0\)、\(C_b[h^V]=-\Box\xi_b=0\)。因此，**不能声称删掉该扇区必然在轨道两侧的光滑真空区域破坏 Lorenz 约束**；\(h-h^V\) 在这两个区域仍可满足约束。只测远离轨道的 divergence 不能验证整个删项操作。

关键在拼接。设 \(k^{ab}=k_+^{ab}\Theta(r-r_p)+k_-^{ab}\Theta(r_p-r)\)，则作为分布

\[
 C^b[k]=C^b[k_+]\Theta_++C^b[k_-]\Theta_-
       +[\bar k^{rb}]\,\delta(r-r_p),\qquad [k]=k_+-k_-.
\]

若完整度规已经满足匹配，逐侧删去 \(h^V\) 后会改变度规的跳跃，并在 \([\bar h_V^{rb}]\ne0\) 时留下

\[
 C^b[h-h^V]=- [\bar h_V^{rb}]\delta(r-r_p).
\]

此时只积分 \((h-h^V):\nabla\nabla\Phi_0\) 会漏掉接触源。具体系数是否非零应由对应模态的跳跃检验；本复核没有把“通常存在”冒充已逐分量测出。即使该特定收缩偶然为零，其他连续性/导数跳跃和 Einstein 源匹配也仍需核验。

这不是仅从一般理论推测：原论文 `:2067–2079` 明说 scalar/vector 用于恢复轨道球面的正则性，并写出分段规范向量产生的额外度规项

\[
 -2[\xi_{(a}]\,\delta^r_{b)}\delta(r-r_p).
\]

**实际消融应保留两种不同解释。** 使用旧的 \(h:\nabla\nabla\Phi_0\) 源并删 vector，可检验“漏项数据经过旧管线，会怎样改变论文可见曲线”；这是历史计算链的指纹试验。它不能自动当作另一个合法、连续、同一粒子源的 Lorenz 解。若解释成规范变换，还必须同时变换 \(\delta\Phi\to\delta\Phi-\xi^a\nabla_a\Phi_0\) 并处理轨道和分布项；仅补上一个 divergence 项也不自动修复被改变的 Einstein 源。

此外，2025 论文 `:626–633` 明确指出其 scalar-only 通量不是该微扰阶的完整规范不变量。故“删了纯规范项，物理通量理应不变”也不能用于这一消融。应比较同一规定规范下的同一定义量。

## 2. 门限错误不是 Kerr 专属；Schwarzschild 是实质反证约束

固定公开 fork 的 `GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl`：

- `:2293–2294` 定义 `lmins1=Max[1,Abs[mm]]`、`lmins2=Max[2,Abs[mm]]`；
- `:2327` 生成 spin −1 径向模式时用了 `lmins2`，`:2422` 拼接该模式时再次用了 `lmins2`；
- `:2437–2438` 的 spin +1 解由前面的 spin −1 数据导出，不能靠这里的正确门限补回已置零的输入。

门限没有 \(a\)，因此在 \(a=0\) 也排除 \(|m_g|=1,\ell=1\)。原 Sam `.m:1209,1375` 使用正确的 `lmins1`；fork 的 `Ret-Numerics/h1-radial-gen` 副本也已使用正确门限（`:2334,2430`）。不能把整个作者仓库等同于单一有错执行路径。

该 vector 结构也不因 \(a=0\) 自动消失：由原论文 `:1314–1318`，取 \(a\omega=0\)、球谐极限 \(b_{\pm1}=I\)，偶宇称角向投影系数变为

\[
 \mathbf S^{(1)}_{l_+l_+}=-2\hat\lambda,
 \qquad \hat\lambda_{\ell=1}=\sqrt2.
\]

这证明没有统一的“乘上 \(a\) 所以 Schwarzschild 无影响”机制；具体激发幅度仍由源匹配决定，不能仅凭角向式宣称所有轨道上的通量影响均相同。

2025 论文补充材料 `main_PRL.tex:681–685` 说其模块化重构数据在 Schwarzschild 极限对所有模态、半径与独立 Lorenz 数据达到机器精度；唯一注明的例外是静态 \((\ell,m_g)=(1,0)\) completion。这里没有说非静态 \(m_g=1\) 使用另一生成器。如果同一有错路径实际用于这项验证，\(m_g=\pm1\) 的数据就应该受到检验。这是需要解释的约束，而不是支持漏项归因的证据。

因此，“Schwarz 对、Kerr 不对”只说明两个实际比较链的差异具有背景/分支/数据依赖；不能直接推出这个无 \(a\) 的门限就是根因。还缺 2025 图对应的代码版本、执行入口及数据谱系。现代公开 HDF5 的缺项指纹，即使很强，也不能补足该历史链条。

## 3. 视界 orbital-energy 因子与本地公式确实一致

令 \(\nu=\omega_c+m_g\Omega_p\)、\(m=m_b+m_g\)、\(k_H=\nu-m\Omega_H\)，\(Z_H\) 为同一物理标量模态在剥去入射相位后的视界幅度，角向范数为 1。则

\[
 N_H=2(r_+^2+a^2)k_H|Z_H|^2
    =4Mr_+k_H|Z_H|^2,
\]
\[
 \dot E^{\Phi,H}=\nu N_H,\quad
 \dot Q^{\Phi,H}_{\rm paper}=-N_H,\quad
 \dot E^{s,H}=\dot E^{\Phi,H}+\omega_c\dot Q^{\Phi,H}
             =m_g\Omega_pN_H.
\]

最后一个等号就是 2025 论文 `:532,542,552`。本地 `environment_cloud.py:100–112` 使用 \(M=1\)、`nh=2*(rp*rp+a*a)*(omega-m*oh)*abs(z_h)**2` 以及 `(omega-omega_c)*nh`，完全对应。其返回的 `charge` 是流入视界的波电荷，符号与论文的“云电荷变化率”相反；没有额外的符号错误。比较的是 effective orbital-energy，不能误用 `wave_energy`，也不能先对各模态取绝对值再相加。

这项解析核验排除了**在相同幅度、角向范数、质量单位约定下**误用 \(\nu\) 代替 \(m_g\Omega_p\) 或漏乘视界面积因子；它没有独立验证 2025 原数据中 \(Z_H\) 的定义和最终单位转换。没有依据修改 \(4Mr_+\) 去拟合曲线。

## 4. 归一化可以排除什么，尚不能排除什么

本地云以直接 Killing 质量积分归一化（`environment_source.py:65–126`），并独立核查 \(E=\omega_cQ\)。现有 `cloud_mass_measures.json` 表明，在保持同一复场系数时，将质量积分测度 \(\Sigma\,dr\,d\Omega\) 改成论文印刷式的 \(r^2dr\,d\Omega\)，\(\alpha=0.3\) 的通量因子仅为 **1.00014896216**。这足以排除“只换这一测度解释约 30% 偏差”，不等于独立验证全部作者归一化。

仍共享/未由原始数值产物独立确认的约定包括：复场应力与电流的共同系数、全角向范数（含 \(2\pi\)）、物理场与 \(\phi^{(1,0)}\) 的 \(\epsilon\) 关系、\(q^2\epsilon^2\) 与 \(q^2M_c/M\) 的图单位转换，以及边界幅度是 \(\Phi\)、\(r\Phi\) 还是其他径向变量。原文 `:454` 的应力式若按通常对称括号权重 \(1/2\) 解读，与 `:515,475–476` 使用的系数存在表面差异；`:462` 的印刷质量式也与本地正 Killing 质量的符号/测度不同。应核对实际实现，不能把排版公式的差异擅自升格成已发现的倍数错误。

若质量和通量始终使用同一个总体复场系数 \(C\)，则固定 \(M_c\) 时 \(\Phi_0\propto C^{-1/2}\)、\(\delta\Phi\propto C^{-1/2}\)，而通量 \(\propto C|\delta\Phi|^2\)，所以 \(C\) 自行消去。现有 `normalization_invariant_flux_ratio.json` 的 \(|F_H|/|F_\infty|\) 在 \(r_p=20M\) 仍与图读数相差 **26.28%**；因此一个共同的云振幅乘数不能单独同时修复两端。该比值仍受有限模态范围和图读误差限制，不能指出究竟哪个源模块出错。

## 5. 本轮消融的新反证：一个模式已足以排除“删项后总通量贴合”

在本文目标的精确阈值 Kerr 云上 \(\omega_c=m_b\Omega_H\)，所以

\[
 k_H=\omega_c+m_g\Omega_p-(m_b+m_g)\Omega_H
     =m_g(\Omega_p-\Omega_H),
\]
\[
 \boxed{\dot E^{s,H}_{\ell m}
 =4Mr_+m_g^2\Omega_p(\Omega_p-\Omega_H)|Z^H_{\ell m}|^2\le0}
 \quad (0<\Omega_p<\Omega_H).
\]

\(m_g=0\) 为零。对相同 \(m\) 的不同 \(\ell\)，实参数角向本征问题的正交性消去交叉项；不同 \(m\) 则由方位积分正交。因此按论文的 scalar-only 定义，其余模式不能以正项抵消一个已过大的负模式。该结论不要求先计算全部模式；也不意味着同一模式内部不同源扇区没有干涉。

已读取同网格消融 `legacy_vector_ablation_alpha03_r20_n88_q18.json`，\(\alpha=0.3,r_p=20M\) 给出

\[
 F^{H,\mathrm{full}}_{00}=-6.7007503411\times10^{-5},\qquad
 F^{H,\mathrm{dropV}}_{00}=-4.8823614232\times10^{-4}.
\]

`figure2_three_orbit_comparison.json` 的同单位图读总量为 \(-5.1763753371\times10^{-5}\)。所以不论删项如何改变其他模式，

\[
 |F^{H,\mathrm{dropV}}_{\mathrm{total}}|
 \ge |F^{H,\mathrm{dropV}}_{00}|
 \simeq9.43\,|F^H_{\mathrm{paper,total}}|.
\]

这在固定其余管线和单位的前提下，已经否定“只漏 vector 就把本地结果推向该论文总通量”。完整 \(00\) 与完整总和不能混称；此处明确利用总和的符号界，而没有把 \(00\) 数值冒充总和。已有独立边界恒等式报告 `legacy_vector_boundary_a03_r20.json` 可交叉核验 vector 振幅；整体误差认证仍由相应数值审计负责。

**原文正 \((2,2)\) 通量的文字范围。** `main_PRL.tex:362` 说 \(r_p\sim50M\) 的 scalar \((2,2)\) 产生正视界通量，却没有在该句明确区分同图中的 Kerr/Schwarzschild。若这句话指上述 threshold Kerr，则与它自己的频率定义 `:471` 和通量式 `:552` 不相容：此时 \(m_g=1\)、\(k_H=\Omega_p-\Omega_H<0\)，连普通 wave-energy \(\nu N_H\) 也为负，不能靠“wave vs orbital”化解。若它指 Schwarzschild，则 \(\Omega_H=0\)、\(\nu>0\)，\((2,2)\) 确可产生正通量，而 \((0,0)\) 的 orbital-energy 可负。故应标为**背景适用范围未讲清、不能套用于 threshold Kerr 的文字**；不能据此认定原代码使用了错误符号，更不能声称它解释现有 32.7% 偏差。

## 建议替换诊断报告的核心结论

“已确认某些公开生成器会漏掉非静态偶极 vector，且现代公开数据带有相应指纹。本地完整度规与原 Sam 方法的抽样对照支持保留该扇区。本轮同定义删项消融使主导视界模态更负约 7.29 倍；阈值 Kerr 的全模态非正恒等式进一步排除了其他模态将其抵消至论文总量。因此只凭此漏项解释 2025 视界偏差的假说不成立，2025 数据谱系仍需确认；删项结果只作为旧计算链消融，不进入生产修正。径向、角向和云形状等已做的精度检查降低了相应数值误差解释的可信度，尚不构成排除所有共享约定、轨道匹配和历史输入差异的证明。”

## 固定来源

- 主文：arXiv:2501.09806v1，`outputs/environment_reference/main_PRL.tex`。
- Lorenz 重构主文：arXiv:2306.16459，`outputs/lorenz_reference/LorenzGaugeKerrCirc.tex`。
- Sam 原实现：`srd24/KerrLorenzCirc` commit `5c1b42793ff893fd0c65a5b48ceecc02d34c5e35`。
- Dyson fork：`ConorDyson/KerrLorenzMSF` commit `6b469f0835e42d3fa531b4997d7a912257119f60`；本文比较有错 `GenerationCodes/Numerics-Asymptotics` 与正确 `Ret-Numerics/h1-radial-gen`，不假定二者都生成过 2025 图。

下列 SHA-256 固定本次实际读取的本地证据；路径中的行号均对应这些字节。

```text
cc1434ceb6eec54b3dbfbf69cbda528bc787ba8e8571eaa9278a2b23bef884b1  outputs/environment_reference/main_PRL.tex
30791a1b1558a5602b04bed43dc44fb79617c345b0e91cdce3a02fba2f44cd05  outputs/lorenz_reference/LorenzGaugeKerrCirc.tex
124fc609fb36f46d4f28378ac040655c8ae7d3cff9346ec62134c358ce223ca9  outputs/paper_metric_reference/srd24_KerrLorenzCirc/metric_reconstruction_calc_radiative.m
0037e4587fe57af9b46d7fc4fc7ad1ba7e521920169c3cfe1ca19d7d762d31a6  outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/GenerationCodes/Numerics-Asymptotics/h1-Radial-gen/MetricReconstructRadiative.wl
2682cf87cf558ff3ada3a0dfa0e8fa3ca12fd73f3f5f069de7d219dccb9c6071  outputs/paper_metric_reference/ConorDyson_KerrLorenzMSF/Ret-Numerics/h1-radial-gen/MetricReconstructRadiative.wl
13f9628f4513c00022ce7ebe69fc3d74db19a8e9fb610ece93ef9095866c8549  src/environment_cloud.py
f8ed45cc7200b9738fa3406d0b8a03552f4847340cabdc92eec97b219c4cebb2  src/environment_source.py
cf398c0cc548f3c88b78d9fcd8333fd1ef49cfa38147fe7d134a846221232324  docs/environment_reproduction/cloud_mass_measures.json
d77d9df592cbae50aeeb62800f7e63eb3a7ef672f7e9ab6e56abf349b5637250  docs/environment_reproduction/normalization_invariant_flux_ratio.json
42aa4d551e8204db9a17a3a4486e990690db336e81ea4ffe1768c0f9cddfc9c2  docs/environment_reproduction/legacy_vector_ablation_alpha03_r20_n88_q18.json
dc959b970fca24293e5ff0a701647cedf4912de8eeb6537aceaafe6cc1ed1370  docs/environment_reproduction/legacy_vector_boundary_a03_r20.json
30ef9b10185960155990ee9bd0213dce5b134aac90c33682f98eabfa3b20c2dc  docs/environment_reproduction/figure2_three_orbit_comparison.json
```
