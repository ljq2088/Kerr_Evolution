# 与参考仓库对齐

基准：[IndigoJX/ParticleSourceEvolutionInHyperboloidalCompactCoord](https://github.com/IndigoJX/ParticleSourceEvolutionInHyperboloidalCompactCoord/tree/13fc2b485f76b1284c2e17879ede9075de062a9d)，提交 `13fc2b485f76b1284c2e17879ede9075de062a9d`。
`code/`、`findings/`、参考 `data/`、`resource_details/`、`README_distribution.md`、`REFERENCE.md` 按该提交原样导入。文件校验值见 `upstream_manifest.json`。这些文件的作者与历史验收声明属于上游；本仓库不另行赋予上游未提供的许可证。

## 坐标与场

新的 s=-2 主实现直接使用上游的 T、R=L²/r、y=-cos(theta)、入射正则方位角 Phi 以及 exp(i m Phi) 约定；内部 M=L=1。演化的是正则四标架中的 peeling 场 psi4，满足 Psi4^R=R psi4。一阶状态是 (P,psi4)，P=A psi4_T+B psi4_R+D psi4。空间为七点有限差分，时间为 classical RK4；轨迹与四块点粒子源也来自上游，完整保留视界、零无穷和轴端点处理。

旧标量实现与 PDF 保留作为 s=0 参考，其 sigma=r+/r、tau=v-2r-4M ln(r/r+) 与上游的关系是

\[
R=\frac{L^2}{r_+}\sigma,\qquad
\tau=T+4M\ln(r_+/M).
\]

后一式要求两边采用同一 tortoise 坐标积分常数；独立初值运行可以分别把初始切片时间记为零。标量场 Phi_scalar=R psi0，所以 psi0=u/L²，u=r Phi_scalar。M=L=1 时这两个振幅数值相同。`src/scalar_upstream.py` 提供坐标和场的显式转换，保存绝对 T 以及从初始切片起算的 T_elapsed，避免混淆时间原点。

## s=0 参考的位置及验证

上游 [RIPLEY_OPERATOR.md](../findings/kerr_point_particle_evolution/theory/appendix/RIPLEY_OPERATOR.md) 给出通用自旋权重算符，可直接取 s=0。本次检查所固定提交的发行文件没有独立的 s=0 演化程序；源项推导中辅助算符的 s=0 标签不能当作标量演化器。

`tests/test_upstream_scalar.py` 独立录入该附录的 s=0 系数，核对时间二阶、混合导数、径向二阶、时间一阶、径向一阶和势项。覆盖 M=1/2.3、a/M=0/0.7/-0.8/0.99、m=0/1/2，共 24 组参数，并包含两个物理端点。两边方程仅差不影响真空解的整体负号。

旧标量求解器仍采用 Chebyshev + 球谐 Galerkin + DOP853。它与上游对齐的是连续方程和坐标/场解释，不能把两者说成相同离散算法；s=-2 路径则直接复用上游算法。

## 复现

```bash
.venv/bin/python -m pip install -r requirements-aligned.txt
bash run_sminus2.sh --run-id kerr-sminus2-smoke --no-progress
.venv/bin/python src/scalar_upstream.py outputs/kerr_m2/evolution.npz outputs/kerr_m2/upstream_coordinates.npz
(cd code/kerr_point_particle_solver && ../../.venv/bin/python -m pytest tests -q)
(cd code/kerr_waveform_tools && ../../.venv/bin/python -m pytest tests -q)
.venv/bin/python -m pytest tests -q
```

默认 s=-2 算例为 a/M=0.8、m=2、质量比 1e-5、128×33 网格，使用上游随仓库给出的通量，无需额外生成 Wolfram/BHPT 通量。它包含点粒子源，与旧的无源高斯标量脉冲是不同物理实验。默认产物为零无穷处的复 psi4_22；并未请求 FFI strain。

## 验证边界

代码与上游相同以及测试通过，不等于完成独立物理验证。保留上游 A1 Weyl 符号/源归一化 bridge 的 `conditional-open` 状态；附录还注明独立 PDF 录入核验未完成。上游取消的高分辨率验证 gates 不记作通过。低分辨率 smoke 演化用于复现和功能检查，不能替代生产分辨率的收敛研究。
