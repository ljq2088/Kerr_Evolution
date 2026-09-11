# 阈值云的视界符号恒等式

依据2501.09806v1补充材料 `eq:scalar_flux`、频率关系和阈值条件，
可直接推导一个不依赖径向求解器的结论。记m=m_b+m_g：

\[
\omega=\omega_c+m_g\Omega_p,\qquad \omega_c=m_b\Omega_H.
\]

于是

\[
\omega-m\Omega_H=m_g(\Omega_p-\Omega_H),
\]

附录的轨道有效标量视界通量化为

\[
\boxed{\dot E^{s,H}_{\ell m}
=4Mr_+ m_g^2\Omega_p(\Omega_p-\Omega_H)|R^H_{\ell m}|^2.}
\]

因此，当0<Omega_p<Omega_H时，每一通道均非正；m_g=0为零。
这个结果对m的正负都成立，与模态振幅、共振强弱、云质量整体
归一化、有限标量模态截断无关。正的整体归一化不改变符号。

另一种验证来自同一公式的能量—角动量关系：

\[
\dot E^{s,H}-\Omega_H\dot L^{s,H}
=4Mr_+(\omega-m\Omega_H)^2|R^H|^2\ge0,
\qquad \dot E^{s,H}=\Omega_p\dot L^{s,H}.
\]

两式同样给出Omega_p<Omega_H时dot E^{s,H}<=0。

v1正文主结果(ii)称在rp约50M附近标量(2,2)给出正的视界通量并
诱导sinking。这与上述附录公式和阈值假设不能同时成立。在本项目
alpha=.3背景，Omega_H约.296293，rp=50的Omega_p约.00282，
确实满足所需不等式。不能将该文字结论直接作为翻转通量符号的依据。

可能需要区分普通波能流、轨道有效能流与未包含的束缚态能量储存/
保守交换；本文档没有声称确定了作者实际代码或图中的具体问题。
这种公式—表述不一致也不自动解释rp20处的32.65%幅值差异。
生产代码继续保留普通波能流、Noether流和轨道有效能流的分别输出。

源文件：`outputs/environment_reference/main_PRL.tex`，标签
`eq:Noether_ind`、`eq:scalar_flux`、`eq:energy_angmom`；原文版本
[arXiv:2501.09806v1](https://arxiv.org/abs/2501.09806v1)。
