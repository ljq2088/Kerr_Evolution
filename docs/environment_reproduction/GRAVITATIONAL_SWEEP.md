# 图4的独立真空引力波基准

固定 a/M=.8771530275949366，计算 rp/M=3.5,4,5 及6到50之间的偶数半径，共26点。每点分别求解 ell=2..12 的全部非零正负 m，共154个通道；m=0静态通道没有辐射。结果单位为每 q²。

每个模式核验 Edot=Omega_p Ldot；汇总时再次核验正负 m 的通量对称性以及直接模态和。每个半径保留BHPToolkit公开的 a=.9、r=10、(2,2)舍入基准检查。全区间从 ell<=10 增至 ell<=12 的最大相对变化约1.04722e-4；这只是有限截断差异，不是严格余项界。

图中同时显示完整有限和与单个正 (2,2) 模式，以避免把单模式曲线当作全部引力辐射。全部原始文件、SHA256、各半径的截断变化保存在 gravitational_flux_sweep_L12.json。此处没有生成尚未齐备的标量/GW比值，也没有消除论文图4已发现的归一化疑问。

重建汇总与图像：运行 src/report_gravitational_sweep.py。底层计算为 src/report_gravitational_flux.py，每个半径显式指定 --radius、--ellmax 12 与 --output。

公开基准：https://bhptoolkit.org/modules/teukolsky/
