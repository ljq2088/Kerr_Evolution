# 辅助质量平方的角向解析导数

这是消除kappa的相近解相减的一项独立准备工作，尚未接入生产Lorenz重构。固定实频率omega，令辅助质量平方为nu。有限球谐基上的标量角矩阵为

\[
H(\nu)=\operatorname{diag}[j(j+1)]-a^2(\omega^2-\nu)C,
\qquad C_{ij}=\langle Y_{im},\cos^2\theta\,Y_{jm}\rangle.
\]

从归一化本征方程H v=A v出发，对nu求导并选取v^T v_dot=0：

\[
\dot A=v^T(a^2C)v,
\qquad
\dot v=\sum_{k\ne L}v_k
\frac{v_k^T(a^2C)v}{A_L-A_k}.
\]

因此S_dot与其theta导数可直接由球谐系数dot v求和。该公式在给定有限基上避免本征解相减；无限基截断误差仍须验证。实现沿用当前的本征向量相位及负m约定，坐标轴处的导数要求单独极限公式，当前入口明确拒绝轴点。

environment_angular_variation.py提供这两个导数入口，生产环境模块尚未调用它。tests/test_angular_variation.py的五项测试通过：与独立四阶中心差分比较，包含ell18/24和负m；验证归一化导数为零、Hellmann–Feynman积分、a=0极限及轴点保护。实际对照误差保存在angular_mass_derivative_audit.json。

接下来还需对径向方程及两端retarded边界求导。若

\[
(\Delta R')'+V(\nu)R=0,\qquad
V(\nu)=K^2/\Delta-\nu r^2-a^2\omega^2+2am\omega-A(\nu),
\]

则

\[
(\Delta\dot R')'+V\dot R=(r^2+\dot A)R.
\]

还必须同时对Wronskian、点源角投影和完整角场求导，才能替换现有trace resolvent的有限差分。单独角向检查不证明径向或完整kappa导数正确，也不证明高阶源匹配问题已解决。
