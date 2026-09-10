# Kerr 点粒子时域演化器

- 来源任务：T007 / W02.3
- 状态：owner accepted
- 公式基准：[Kerr 点粒子 Teukolsky 演化](../../findings/kerr_point_particle_evolution/README.md)
- 数值合同：[numerical_contract.md](docs/numerical_contract.md)
- 正式数据：[SXS:BBH:0305 medium bundle](../../data/kerr_point_particle_evolution/sxs0305_medium/README.md)

该包演化固定 $m$ 的双曲紧致化 Teukolsky 方程，包含 Ori--Thorne 初值、唯一
$(R_p(T),\Phi_p(T))$ 轨迹、七点空间离散、四块点粒子源、classical RK4、
$\mathscr I^+$ 投影、checkpoint/restart 和 FFI strain。Kerr 背景、轨迹、波形 I/O
与 strain 的公共实现位于 [`kerr_waveform_tools`](../kerr_waveform_tools/README.md)。

## 正式运行

仓库内不显式给出输出目录时，`run_case.py` 默认写入
`data/kerr_point_particle_evolution/runs/.staging/<run_id>/`；全部请求的数据、图、
hash 和 loader 检查通过后，目录原子定稿为 `runs/<run_id>/`。初值按 versioned key
缓存于 `data/kerr_point_particle_evolution/initial_data/`，同一物理参数再次运行时复用。

~~~bash
PYTHONPATH=code/kerr_point_particle_solver/src:code/kerr_waveform_tools/src \
python code/kerr_point_particle_solver/scripts/run_case.py \
  code/kerr_point_particle_solver/configs/smoke/evolution.json \
  code/kerr_point_particle_solver/configs/smoke/output.json \
  code/kerr_point_particle_solver/tests/fixtures/chi0p8_l6_flux.json
~~~

Production output 默认应同时请求 complex $\psi_{4,\ell m}$、$H_{\ell m}=rh_{\ell m}$
及两者的实部/绝对值图。`status=complete` 表示计算和交付完整，owner 审阅状态由
bundle 中独立字段记录。

## 检查

~~~bash
PYTHONPATH=code/kerr_point_particle_solver/src:code/kerr_waveform_tools/src \
python -m pytest code/kerr_point_particle_solver/tests
~~~

最低分辨率跨代码对比已由 owner 接受；原计划中其余数值 gates 被取消而非记为通过。
FFI 的 A1 sign/normalization bridge 继续在 artifact metadata 中标为
`conditional-open`。
