# 静态圆轨道局部匹配：实现与未通过项

本项仍是诊断，不是完整静态 Lorenz 点粒子解。物理边界尚未施加，
Kerr 的全部局部跳跃条件尚未收敛，不能接入论文尾迹或总通量验收。

## 圆轨道对称性修正

单支 IRG 经 Lorenz 变换并加复共轭后，虽然具有正确的两种极端 Weyl
曲率，仍可包含违反圆轨道对称性的规范部分。Kerr 度规及静态圆轨道
应力张量在 I:(t,phi)->(-t,-phi) 下不变，因此采用

    h_circ = (h + I* h)/2,
    (I* h)_ab = p_a p_b h_ab, p=(-1,1,1,-1).

这是张量拉回的平均，不是逐点削减误差。它消去 tr、t theta、r phi、
theta phi 分量。线性化 Einstein 算子和 Lorenz 条件与背景等距变换
可交换，故真空解及规范条件得以保持。独立曲率计算还验证圆轨道
源的 Psi0/Psi4 保持不变：在 a=0/.6、ell=2/3、轨道内外的测试点，
相对差约 1e-15。未对未验证的一般源作此假设。

实现位于 sourced_static_spin2；circular_symmetry=False 保留原始分支
供对照。静态 spin-2 底层函数仍返回原始复 IRG/Lorenz 分支。

## 匹配方程

particular 度规是 ell>=2 的 sourced spin-2 与偶数 ell trace 之和。
各分支 chi、kappa、B 的 particular 径向数据仍在 r0 取零。
加入自由标量 Hessian 及 B,C,D,E,F,G 补全：

    h_free,ab = nabla_a nabla_b kappa,
    kappa = R_j(r) Y_j0(theta),
    (Delta R_j')' - j(j+1)R_j = 0.

每个偶数 j 使用 r0 处 (R,R')=(1,0),(0,1) 两个局部基。
由于 Ricci=0 及 Box kappa=0，该 Hessian 无迹且满足真空 Lorenz 方程。
系数只表示 outside-minus-inside 跳跃，不能解释为最终 Up 振幅。

所有十个协变分量采用光滑角向向量 sin(theta) partial_theta 投影，
测试函数是 P_j(cos theta)。两侧 order-8 Jets 在 r0±epsilon 求值，
通过二阶 Taylor 继续到轨道。目标是

    [h_ab]_j=0,
    [partial_r h_ab]_j = -8 P_j(0)/(ut Delta0)
                            (u_a u_b + g_ab/2).

常规诊断同时拟合这两个条件；它不能独立验证导数跳跃。
--continuity-only 则只拟合场连续性及 Delta c_E=E、Delta c_G=L-aE，
把所有导数跳跃保留作独立验证。它目前以最小二乘加入守恒量行，
并非代数消元强制精确成立，因此输出同时保留拟合值和期望值。

报告缩放为 h_ab/(L_a L_b)、r0 partial_r h_ab/(L_a L_b)，
L=(1,1,r0,r0)。这是坐标分量诊断尺度，绝不是通量误差界。
列范数归一化前先排除几乎没有投影的列，避免把舍入噪声放大。

## 结果与限制

全部在 r0=6、epsilon=5e-5。

* a=0、L=6、12点角积分，拟合 j<=4：缩放残差 8.87e-10。
  E/G 未施加而独立恢复至约 4e-10/2e-10 绝对差。
  该小残差包含用于拟合的导数条件，不能当作独立导数验证。
* a=.6、L=6、free j<=6、12点、仅连续性拟合 j<=4：
  低阶 j<=2 的值/独立导数残差分别 2.67e-6 / 6.66e-6。
* a=.6、L=8、free j<=8、16点、仅连续性拟合 j<=6：
  同样 j<=2 的值/独立导数残差为 7.09e-8 / 9.02e-7；
  j<=4 为 1.97e-6 / 4.34e-6。
  但最高拟合阶附近的导数残差仍达 5.90e-2，因此整体未通过。

这两次 Kerr 试验同时改变截断、角积分和拟合阶数，只能说明低阶
结果改善，不能独立量化其中任一误差来源。还要固定物理投影阶数
提高各项分辨率，检查补全系数、自由标量自由度和 Taylor 外推。
让过高的自由标量阶数只受低阶投影约束可出现显著病态；对应失败
试算也保存，不能从小的拟合连续性残差推断物理解正确。

原始未施加圆轨道对称性的 L4 试算标记为 historical，保留它说明
为什么真空验证与 Weyl 验证本身不足以证明完整点粒子匹配。

即使局部条件全部通过，还必须把 particular 的全局齐次自由度
适配到视界/无穷远，并明确采用正确质量角动量解还是 Berndtson 解。
参考半径为零数据的当前补全基不能直接套用论文的边界适配系数。

复现入口：src/report_static_matching.py。数据为同目录 static_matching_*.json。
新增独立对称性/自由标量测试通过，主项目 tests 共 62 passed。
