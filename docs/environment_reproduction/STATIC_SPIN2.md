# 静态spin-2 Hertz–Lorenz构造

实现依据2306.16459v3的零频构造，使用m=omega=0专用方程，不取含
1/omega的非静态表达式的数值极限。当前已得到真空分支及点粒子曲率
归一化，但还没有完成静态各扇区的联合匹配和物理边界适配。

## 径向和角向基

取lambda=ell(ell+1)、x=(r-1)/sqrt(1-a²)，R0选Legendre P_ell(x)
或Q_ell(x)，满足(Delta R0')'-lambda R0=0。静态Hertz径向函数为

    P=Delta² R0'',
    P'=(lambda-2)Delta R0',
    P''=lambda(lambda-2)R0,
    Delta P''-Delta'P'-(lambda-2)P=0.

角函数为单位球面范数的spin=-2、m=0球谐：

    Y_-2,l0 = sqrt[(2ell+1)/(4pi) * (ell-2)!/(ell+2)!]
              * sin²(theta) * d²P_ell(cos theta)/d(cos theta)².

psi=P Y_-2,l0。通过文献IRG算符构造复度规，再加上其复共轭形成
实静态度规。单独的复Hertz构造带有手征性，不能把其某一组为零的
extreme Weyl分量误判为整个物理扰动纯规范。

## 零频规范向量与符号核对

定义rho=r+i a cos(theta)、rhoc=conjugate(rho)，
U^{ab}=l_+^a m_+^b-m_+^a l_+^b，L_n=partial_theta+n cot(theta)。
在本实现的显式约定下，使用

    rhoc H_U = -[rhoc partial_r L_2 -2(L_2+i a sin(theta)partial_r)]psi/(2lambda),
    H^{ab}=H_U U^{ab}/Sigma,
    xi^a=rhoc² nabla_b H^{ab}-nabla^a chi,
    h_L,ab=h_IRG,ab-2 nabla_(a xi_b).

H_U前面的负号与原文所印零频表达式不同，但由它前面的有源spin-1
方程固定，而非用残差拟合。例如a=0、ell=2时，R0=P2(r-1)、
P=3Delta²、lambda=6。若取原文所印正号，令

    F=(rP'-2P)/(2lambda)=r³(r-2)/2,

直接求导得到(Delta partial_r²-lambda)F=+P，而H_U的有源方程要求
-P。取负号后方程成立，且完整度规的Lorenz条件与无迹条件同时恢复。
该明确代入例子保留此符号差异，便于与作者约定继续核对。

chi满足(Delta chi')'+angular_Laplacian chi=S。将Y=Y_-2,l0记为
背景角因子，源为

    Re S=-(rP'-2P)L_1 L_2 Y/lambda
          +a² P''[sin(theta)cos(theta)L_2Y+2sin²(theta)Y]/lambda,
    Im S=a[P'cos(theta)L_1L_2Y+rP''sin(theta)L_2Y]/lambda.

实部仅耦合j=ell-2,ell,ell+2，虚部仅耦合j=ell-1,ell+1。各项投影
到标量Y_j0后，解(Delta chi_j')'-j(j+1)chi_j=S_j。当前在参考半径
指定chi_j=chi_j'=0，保留待匹配的齐次规范自由度。

## 点粒子曲率归一化

静态physical实度规的曲率关系为

    Psi0 = (1/2) partial_r^4 conjugate(psi_Hertz).

利用P''''=lambda(lambda-2)P/Delta²，把独立pybhpt静态点粒子Psi0
的In/Up径向解与Legendre分支相匹配，得到Hertz振幅A。物理度规使用
A h_L+conjugate(A h_L)。归一化在一个真空半径完成，随后在别的半径
独立检查Psi0和Psi4，两者均恢复输入点粒子曲率。

测试覆盖Schwarzschild/Kerr、P/Q分支、ell=2/3及轨道内外。真空方程、
Lorenz条件和迹按度规幅度缩放后检查；另检查实IRG与实Lorenz度规的
extreme Weyl分量不变且非零。点粒子两种曲率的独立比较相对差约1e-15。
全套主项目测试60 passed。这些检查不证明跨轨道连续性或最终边界条件。

下一步必须联合static trace、补全基及自由标量齐次项，确定轨道处的
跳跃系数并验证所有度规分量。不能直接将当前单独spin-2分支当作完整
环境源，也不能直接采用未适配到本基的补全系数。

后续联合匹配发现单支IRG还保留圆轨道对称性为奇的规范部分。
sourced_static_spin2现默认以(t,phi)->(-t,-phi)等距拉回作平均，
独立确认两种极端Weyl曲率不变。修正与尚未通过的Kerr局部条件
详见 `STATIC_MATCHING.md`；底层复分支及其推导保持原样。
