# Kerr waveform 公共工具

- 来源任务：T007 / W02.3
- 状态：owner accepted
- 物理基准：[Kerr 点粒子 Teukolsky 演化](../../findings/kerr_point_particle_evolution/README.md)

该包保存 W02 solver 与 W03--W06 可以共同使用的实现，不依赖完整时域演化器：

- `background.py`：Kerr 几何、ISCO、light ring、tortoise 与紧致化坐标；
- `transition.py`：Ori--Thorne 输入、初值生成和 versioned cache；
- `worldline.py`：$T$ 参数化轨迹、事件与 stage-time jets；
- `harmonics.py`：spin-weighted spherical harmonic 与 mode 投影；
- `waveform_io.py`：complex $\psi_4$ 的 v1/v2 兼容读取；
- `strain.py`：完整时间序列上的 FFI 和 $H_{\ell m}$ v1/v2 兼容读取；
- `artifacts.py`、`plotting.py`：原子 artifact、hash 和基础波形图。

`scripts/` 提供 ISCO flux request、初值 cache、波形/strain 绘图和 manifest 检查。
包内不包含 fixed-$m$ PDE、source、网格、RK4 或 checkpoint。

~~~bash
PYTHONPATH=code/kerr_waveform_tools/src \
python -m pytest code/kerr_waveform_tools/tests
~~~
