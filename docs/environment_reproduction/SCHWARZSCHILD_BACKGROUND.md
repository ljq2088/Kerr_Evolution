# Schwarzschild比较所需的准束缚背景

已实现复频率|211>背景的独立双向径向shooting，并构造复径向函数及
协变Hessian。尚未接入Schwarzschild环境通量求和，不构成图2/3验收。

对于a=0、ell=1，径向方程为

    (Delta R')' + [r^4 omega²/Delta - mu² r² - 2] R = 0,
    Delta=r(r-2).

视界使用exp(-i omega r*)乘四阶Frobenius展开；远区使用
exp(ikr)r^p乘六阶逆r展开，其中Im(k)>0，选择空间衰减支。
在r_match=2/mu²匹配两侧对数导数，同时求实部和虚部为零。
在通用二维求根后检查残差，必要时以复Newton修正小的阻尼率，
避免仅由大得多的频率实部决定停止条件。

| alpha | Re(M omega) | Im(M omega) |
|---|---:|---:|
| 0.2 | 0.198952661253 | -4.06043825e-8 |
| 0.3 | 0.296192346399 | -9.45565167e-6 |

将远边界从45个氢原子近似衰减长度推到60个，并将视界距离1e-5降至
1e-6后，显示精度内的频率相同，匹配残差约1e-13或更小。这里的
零频率差不等于严格误差为零；两次计算都使用同一双精度shooting方法。

## 在穿越视界的切片上定义云质量

衰减的准束缚态不能直接照搬静态Kerr阈值云的BL定时切片归一化。
使用ingoing Kerr-Schild时间T=t+h(r)，h=2 log((r-2)/2)。
令R_KS=exp(i omega h)R_BL，视界处的对数相位相消，场正则。
其径向导数为

    R_KS'=exp(i omega h)[R_BL'+i omega h' R_BL].

逆度规非零时间径向分量为g^TT=-(1+2/r)、g^Tr=2/r、g^rr=1-2/r。
角函数取单位球面范数后，直接应力张量和Noether流积分给出

    E = integral r² [(1+2/r)|omega|²|R_KS|²
                      +(1-2/r)|R_KS'|²
                      +(ell(ell+1)/r²+mu²)|R_KS|²] dr,
    Q = integral 2r² [(1+2/r)Re(omega)|R_KS|²
                      +(2/r)Im(R_KS* R_KS')] dr.

背景按T=0的E=1归一化。相应吸收率为

    F_H^E=2r_H²|omega|²|R_KS(r_H)|²,
    F_H^Q=2r_H²Re(omega)|R_KS(r_H)|².

未使用这些平衡式进行归一化，而将其作为独立检查：

    2 Im(omega) E + F_H^E = 0,
    2 Im(omega) Q + F_H^Q = 0.

当前以r_H+1e-6近似读取视界振幅，alpha=.2/.3两组的相对平衡残差
均小于2e-6。精确复频率背景在多个外部半径的协变KG残差通过检验。
这种有限偏移的平衡验证不能称为无误差的视界极限。

## 把频率虚部置零的明确含义

可选freeze_decay=True保留复本征频率求出的径向函数，只将时间因子
替换成exp(-i Re(omega)t)。在BL坐标中，冻结后的场满足

    (Box-mu²)Phi_frozen
      = [Re(omega)²-omega²]/(1-2/r) * Phi_frozen,

而不是精确齐次KG方程。该非零残差已通过直接Hessian收缩核对。
质量归一化仍指底层复频率物理解在KS T=0切片上的能量，不能将
冻结后的近似场宣称为严格稳态、严格守恒的云。

这是对“固定径向本征态而冻结时间衰减”的明确实现。论文说明置零虚部
导致背景不再精确满足KG方程，但具体数值归一化/冻结操作仍须核对。
不能在完成核对前宣称本实现复现了作者的Schwarzschild通量。

可复算入口：src/environment_schwarzschild_cloud.py。
证据：schwarzschild_cloud_spectrum.json、schwarzschild_cloud_balance.json、
tests/test_schwarzschild_cloud.py。

本轮主项目测试57 passed；其中两项新增测试分别覆盖精确复频率KG与
能量/荷衰减平衡，以及冻结时间频率后非零KG缺陷的解析表达式。

## 受迫通道接入

采样器支持--background schwarzschild-frozen，记录复本征频率、KS切片
质量归一化及冻结时间衰减的近似。它沿用论文的h:H源算符，不把背景
宣称为精确稳态，也不宣称完整通量严格守恒。

alpha=.3、rp=20、标量ell=m=3，径向每面板4点、角向10点、
ell_g<=4、首面板16个log点的首轮结果为：无穷远有效通量
1.457479072e-5，视界有效通量8.541997910e-13（剥离q² eta）。
这是单通道，不能与论文总通量直接比较。径向8点加密正在计算。

球对称性提供独立检验：背景ell_b=1时，标量ell=3只允许
ell_g=2,3,4参与；标量ell=4,m=3被赤道奇偶规则禁止。
在r=10逐ell投影，超出三角条件的相对贡献为7.16e-12，禁戒奇偶
贡献为7.28e-11，见schwarzschild_source_selection.json。它仍是单点诊断。
