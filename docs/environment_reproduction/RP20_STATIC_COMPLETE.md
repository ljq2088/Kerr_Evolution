# rp=20M的五个静态标量通道与完整有限场

alpha=.3、a=.8771530275949366、Lg=18的静态度规源采样已经完成。
最初的1000M Green边界批次成功计算ell=3,5,7,9；ell=11的逆半径
展开末项相对量.00218063698超过1e-3保护阈值，因此该批次失败，
四个已完成通道保留。它们不是被设为零的缺失场。

把标量Green外边界移到4000M后，ell11末项指标为1.50880595e-6，
五个静态通道均完成。再将边界移至8000M，逐项确认所有88个源
节点、权重和复源值相同；五个通道在粒子位置的场值相对变化
依次为9.17e-12、2.18e-11、1.97e-10、9.65e-12、4.88e-12。
公共径向采样点最大场值变化除以最大原场值也均小于2e-10。
阈值静态通道的轨道有效通量严格为零，但场值非零。

比较的是物理径向场和视界振幅，不比较束缚Up解的任意归一化。
该检查只涉及相同源下的外边界移动，不是源网格、静态匹配、
度规截断或整个连续区域的收敛证明。

## 版本一致的恢复

原始度规缓存的数值源码哈希为
45328186c7c344105d71d58fde1df86eba6e2d0710d4b0667e12a52e42f9fddb，
对应Git提交66da5c7e1239be84ae751d8462c7bb1e146f18cd。当前源码已有
后续缓存优化，严格缓存键不同。最初启动的当前版本4000M重算
已主动停止；没有把历史缓存改写为当前版本。

通过report_archived_mode_batch.py从上述提交提取数值源码，先
验证源码哈希，再使用历史版本模块读取其自己的缓存。仅输出
路径指向当前项目。4000M与8000M两个批次均记录88个缓存命中、
零个重新构造的径向度规点；缓存元数据和张量值未被修改。
版本记录为static_r20_recovery_source.json。

完整批次：scalar_batch_mg0_1b048d12df97.json（4000M），
scalar_batch_mg0_677d2a53818c.json（8000M）。逐通道比较及输入
哈希见static_r20_boundary_4000_8000.json。

8000M检查的可重复命令（在仓库根目录运行；4000M仅改边界参数）：

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python src/report_archived_mode_batch.py \
  --source-commit 66da5c7e1239be84ae751d8462c7bb1e146f18cd \
  --expected-source-hash 45328186c7c344105d71d58fde1df86eba6e2d0710d4b0667e12a52e42f9fddb -- \
  --metric-m 0 --scalar-ells 3 5 7 9 11 --metric-ellmax 18 \
  --radial-order 8 --angular-order 18 --horizon-log --horizon-order 32 \
  --source-inner-offset .0005 --source-outer-radius 320 --workers 2 \
  --static-matching docs/environment_reproduction/static_tetrad_a0.877153_r20_L18_q32_j18_free18_paper_eps5e-06.json \
  --green-outer-radius 8000
```

## 当前结果与未完成项

使用4000M静态响应的新清单为
flux_coverage_L18_nt18_i12_h12_f12_m0go4000_m1go4000_m2go4000.json。
场ell=2..12达到88/88通道；无穷远36/36及非零视界85/85有限通量
总和保持原值，因为补入的静态通道不贡献轨道有效通量。

particle_field_L18.json现在包含全部11个角阶的完整m和。
新增奇数阶的局部尖点诊断与偶数阶平滑衔接：
-L^2 B_ell/(epsilon q)在ell=9、11的实部分别为.3219824069、
.3132248896，独立局部预测为.2976938312。有限阶趋近不等于
大ell收敛证明，且该诊断仍区分球面和椭球面模态。

完整场只属于rp=20M的有限分辨率结果；论文图1要求rp=3.5M，
其队列仍在运行。图7的幅度差异、静态匹配截断边缘、全源收敛
及论文其余半径/参数的结果尚未解决。
