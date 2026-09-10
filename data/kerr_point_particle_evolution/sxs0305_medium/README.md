# SXS:BBH:0305 medium 点粒子波形

- 状态：`complete`；`review_status=owner_accepted`
- 权威清单：[manifest.json](manifest.json)
- 背景：$M_f/M=0.952032939704$，$\chi_f=0.6920851868180025$
- 网格：$(N_R,N_y)=(512,129)$；$\Delta T_{\rm out}=0.1M$

`mode_22/` 与 `mode_44/` 各保存完整 complex $\psi_{4,\ell m}$、FFI
$H_{\ell m}=rh_{\ell m}$、实部/绝对值图和原始 run metadata。`input/` 保存 SXS
$h_{22},h_{44}$、remnant metadata、ISCO flux、plunge 初值和两个演化配置。

完整性以 `manifest.json` 中的逐文件 SHA256 为准，可用
[`validate_file_manifest.py`](../../../code/kerr_waveform_tools/scripts/validate_file_manifest.py)
复核。原始 $\psi_4$ 是权威演化结果；strain 是保留 provenance 和
`conditional-open` A1 convention 状态的派生产品。
