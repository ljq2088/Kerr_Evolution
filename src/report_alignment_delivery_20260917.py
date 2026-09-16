"""Assemble measured alignment status without changing physical results."""
from pathlib import Path
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'docs/field_alignment_20260917'
def read(name):return json.loads((OUT/name).read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 field=read('li_field_alignment.json');flux=read('li_new_radius_flux.json')
 bc=read('scalar22_outer_boundary_control.json');li4=read('scalar22_Li_order4_boundary_control.json')
 rows=[]
 for r in field['ring_comparisons']:
  rows.append(dict(orbit=r['orbit'],mode_selection=r['selection'],observer_radius=r['radius'],angular_rms_ratio=r['local_angular_rms']/r['reference_nominal_angular_rms'],normalized_amplitude_rms_difference=r['normalized_amplitude_rms_difference'],fraction_within_reference_color_interval=r['inside_color_interval_fraction']))
 with (OUT/'alignment_ring_summary.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 chosen=[('coulomb',None,32000),('coulomb',None,64000),('Li_exact_prefactor_series4',4,4000),('Li_exact_prefactor_series4',4,8000),('Li_exact_prefactor_series4',4,32000)]
 labels=['Coulomb\n32,000','Coulomb\n64,000','Li form, order 4\n4,000','Li form, order 4\n8,000','Li form, order 4\n32,000']
 by={(r['orbit'],r['boundary_method'],r['inverse_r_order'],r['outer_radius']):r for r in bc['rows']+li4['rows']}
 fig,ax=plt.subplots(1,2,figsize=(12,4.7),layout='constrained')
 for radius,marker in [(50,'o'),(150,'s')]:
  values=[100*next(v['normalized_amplitude_rms_difference'] for v in by[(42.1,*c)]['full18_only_scalar22_replaced'] if v['radius']==radius) for c in chosen]
  ax[0].plot(range(len(chosen)),values,marker=marker,label=f'r = {radius} M')
 ax[0].set_title('Bound branch: orbit 42.1 M');ax[0].set_ylabel('Amplitude RMS difference from Li image (%)')
 for key,label,marker in [('nominal_total_over_Li','Nominal asymptotic flux','o'),('actual_current_total_over_Li','Conserved radial current','s')]:
  ax[1].plot(range(len(chosen)),[by[(41.1,*c)]['total_infinity_flux_only_scalar22_replaced'][key] for c in chosen],marker=marker,label=label)
 ax[1].axhline(1,color='black',lw=.8,ls='--');ax[1].set_title('Radiative branch: orbit 41.1 M');ax[1].set_ylabel('Total infinity flux / Li plot value')
 for a in ax:
  a.set_xticks(range(len(chosen)),labels,fontsize=9);a.grid(alpha=.2);a.legend(fontsize=9);a.axvspan(1.5,len(chosen)-.5,color='orange',alpha=.07)
 fig.suptitle('Published fourth-order boundary: fixed-source cutoff convergence control',fontsize=12)
 fig.savefig(OUT/'boundary_control_summary.png',dpi=190);fig.savefig(OUT/'boundary_control_summary.svg');plt.close(fig)
 report=r"""# 明显差异对齐核验（2026-09-17）

本轮重新计算了两个匹配轨道半径的 36 个响应通道，并统一了场变量、模态集合、参考图颜色映射和取样坐标。辐射支的幅度尺度已与后续独立结果接近；束缚支仍有显著节点差异，不能宣布完整复现，也不能把所有剩余差异归为普通数值误差。

## 计算与比较范围

参考为 [Li 等，arXiv:2507.02045v2](https://arxiv.org/html/2507.02045v2) 的标量场图及通量图。该文明确指出 Dyson 原文的标量分离常数和视界通量归一化需要修正。它是后续独立参考，不是 Dyson 原图的同一组原始数值。

本轮采用 α=0.3，轨道半径 rp/M=41.1、42.1，观察切面为 t=0、θ=π/2。场量为 |δΦ|/(q α³√η)，η=Mc/M。主结果包含 2≤ℓ≤5 的全部 18 个允许模式；另保留论文正文“贡献无穷远通量模式”的六个正 m 模式解释。两组均未调幅或调相。

局部背景是精确同步云 a/M=0.8771530275949366，参考标注 a/M=0.88。度规截断 L≤6，角向阶数 12，径向阶数 8，近视界对数分辨率 32。本轮是有限分辨率核验，不是已收敛的同背景精密复现。全部 36 个源通道由当前算法重算，带来源哈希。

## 已纠正的对照问题

1. 不能将原图色条上限当作场的全域最大值；原图中心含有色条外颜色。现在使用原始颜色查找表，在相同坐标取样。
2. 不再把 rp=41.6/41.8 的旧计算直接与后续 rp=41.1/42.1 的图比较。
3. 原始 Dyson 图 5 在固定位置仍有约 3.7–4.4 倍差异；图 1 的偏差方向相反，因此一个全局归一化因子无法同时解释。
4. 当前场峰值为 0.362366、0.371678，与 Li 图色条的 0.37 同量级。这只是幅度尺度一致的证据，不能以峰值代替逐点检验。

## 逐圈定量结果

定义 D(r)=sqrt(mean[(|Φlocal|−|Φref|)²])/sqrt(mean[|Φref|²])，在同一圆周的 72 个角度计算。参考值来自图像色带；色带区间不是统计误差，也不是作者原始数组。

| 轨道 rp/M | 观察半径 r/M | 场振幅 RMS 差异 D | 判断 |
|---:|---:|---:|---|
"""
 for rp in [41.1,42.1]:
  for radius in [50,70,100,150,180]:
   r=next(x for x in rows if x['orbit']==rp and x['mode_selection']=='All 18 allowed modes' and x['observer_radius']==radius)
   report+=f"| {rp} | {radius} | {100*r['normalized_amplitude_rms_difference']:.2f}% | {'明显偏差' if r['normalized_amplitude_rms_difference']>.3 else '仍有残差'} |\n"
 report+='\n无穷远总轨道能通量（相同标量 ℓ≤5，单位 q²η）：\n\n| rp/M | 本地 | Li 图读数 | 相对差异 |\n|---:|---:|---:|---:|\n'
 for r in flux['rows']:
  v=r['references']['Li'];report+=f"| {r['orbit']} | {r['local_infinity_flux']:.9e} | {v['reference']:.9e} | {v['relative_difference_percent']:+.2f}% |\n"
 report+=r"""
41.1 处 +21.82% 的差异是本轮新发现的未解决项，不能因远离阈值的通量较吻合而略去。参考通量来自矢量图折线插值，仍不等同于作者原始数据。这里没有计算完整视界总通量，因为这些场图模式缺少 ℓ=0、1。

## 最有辨别力的算法控制

固定源 J 和其余 17 个模式，只改变临界附近标量 (2,2) 的外边界算法。先使用 Li 正文的精确前因子、独立推导的四阶展开：

| 外半径 Rb/M | 42.1 的 r=50 场 RMS 差 | 42.1 的 r=150 场 RMS 差 | 41.1 名义总通量/Li | 41.1 真实守恒总流/Li |
|---:|---:|---:|---:|---:|
| 4000 | 20.67% | 14.97% | 0.7471 | 1.1995 |
| 8000 | 63.14% | 59.53% | 1.1547 | 1.2189 |
| 32000 | 63.14% | 59.53% | 1.2181 | 1.2182 |

四阶 4000M 的场更接近参考，但它不收敛；移到 32000M 后，场与库仑方法相符至约 3×10⁻⁶，名义总通量相差约 0.012%。这些是预先指定的截断点，没有搜索一个“最佳贴图半径”。作者公开 notebook 的 B1–B4 已用独立解析器在 80 位精度代入，与独立递推共八项的最大相对差为 1.43×10⁻⁴⁵；系数一致性已闭合。作者仓库固定版本为 `da78d2c96af5854e0feea0421af1ebb6c45c6bc0`，仅含边界系数，没有求解驱动器及原图端点设置。

辅助零阶/六阶控制：

- 库仑/Whittaker 边界从 32,000M 移到 64,000M，两个轨道的取样场变化小于 3×10⁻⁹。
- 故意只用零阶渐近边界、放在 4,000M：42.1 的 r=50、150 圆周差异分别从 63.14%、59.53% 降至约 13.76%、13.64%。移到 8,000M 后又返回已收敛基准。表面接近参考并不代表正确边界。
- 同一个零阶 4,000M 设置，在 41.1 处按“单位无穷远振幅”公式得到的总通量只比参考高 3.05%；但 (2,2) 模的实际守恒流是该名义通量的 1.4939 倍，总守恒流仍比参考高 19.31%。这揭示了边界幅度归一化失真。

这里控制展开的参数是 k=√(ω²−μ²)，并非只看 ω 是否较大。作者公开 B1 在 k→0 时含 k⁻³ 项；按本地两个算例，在 4000M 处 |B1/R|≈3.41、4.53，首修正已大于首项。

因此，近阈值长程库仑区内过早施加渐近条件，确实能同时制造“节点偏移”和“通量看似对齐”。这是明确的算法敏感性及一条具体根因候选；Li 正文明确给出四阶展开，并公开了边界系数 notebook；未给出生成原图的外截断半径及原始数组。上述零阶/六阶试验并非作者四阶算法的精确复刻，目前不能断言他们实际采用了此设置。不能把这种未收敛解作为修复后的物理解。

## 已排除或限制的解释

- 角向函数经独立球谐、SciPy 椭球函数和投影检验；未发现相位约定错误。
- 独立 Gauss 积分与累计 Green 构造对 R、R′ 的差异小于 1.6×10⁻⁹，最终相位合成差异小于 3×10⁻¹⁰。
- 旧 rp=41.6 的当前源三点复算与缓存相符；这是有限抽查，不能代替全域收敛。
- 近视界遗漏壳层的控制仅改变主模视界通量约 0.134%，解释不了场的倍数差异。
- 固定 a=0.88 的云频率由复 Leaver 与独立双端射击一致算得 Mωc=0.29629353472561146+2.2166094×10⁻⁹i。把精确实部与 a=0.88 放入 Green 传播器后，束缚支 r=50 的差异约从 63% 降到 58%，仍未解决。该试验冻结原来的度规源及云剖面，不是背景一致的新计算。
- Li 源函数图的径向形状接近（约 1.4% 内），但图示绝对尺度约 21.6 的常因子尚未从作者作图约定中闭合；不能据此擅自缩放物理源。

## 交付内容与复算

- `li_field_alignment.png`：原图与两种模态范围的本地场。
- `li_field_ring_comparison.png`：相同圆周上的角向分布。
- `boundary_control_summary.png`：外边界控制及真实守恒流对照。
- `alignment_ring_summary.csv`、`li_field_alignment_points.csv`：定量取样表。
- `li_field_alignment.json/npz`：场数组、统计量、来源哈希。
- `nonstatic/` 与四个 `forced_mode...mg0...json`：本轮 36 个源响应通道。
- `MATHEMATICAL_CONVENTIONS.md`：方程、归一化和源/Green 对应关系。

在仓库根目录设置 `OPENBLAS_NUM_THREADS=1`、`PYTHONPATH=src`。非静态计算入口为 `src/report_li_field_alignment_runs.py`，静态入口为 `src/report_li_static_new_orbits.py`。绘图入口 `src/report_li_field_alignment_plot.py` 的 `--extra` 应传入本轮四个 mg0 静态响应文件；完整输入和文件哈希见结果 JSON。边界控制入口为 `src/report_li_scalar22_outer_boundary_control.py`；本报告入口为 `src/report_alignment_delivery_20260917.py`。

本轮修正了对照流程并补齐新算例，没有给物理求解器加入拟合因子或未经证明的公式替换。仍待闭合的是：背景一致的 a=0.88 准束缚云源、完整分辨率控制，以及参考计算中未公开的外边界与源图归一化细节。当前不能宣称所有明显差异均已消除。
"""
 (OUT/'README.md').write_text(report,encoding='utf-8')
 validation={'status':'partial_alignment_with_explicit_unresolved_discrepancies','production_physics_modified':False,'fresh_response_channels':len([p for p in field['inputs_sha256'] if '/nonstatic/' in p or '/forced_mode_' in p]),'amplitude_fitted':field['amplitude_fitted'],'phase_fitted':field['phase_fitted'],'inputs_sha256':{n:digest(OUT/n) for n in ['li_field_alignment.json','li_new_radius_flux.json','scalar22_outer_boundary_control.json','scalar22_Li_order4_boundary_control.json','fixed_source_exact_spin088.json']},'implementation_sha256':digest(Path(__file__))}
 assert validation['fresh_response_channels']==36
 for name,expected in field['inputs_sha256'].items():
  assert digest(ROOT/name)==expected, name
 for name,expected in field['assembly_dependencies_sha256'].items():
  assert digest(ROOT/'src'/f'{name}.py')==expected,name
 validation['field_input_and_dependency_hashes_verified']=True
 (OUT/'delivery_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
 print(json.dumps(validation,indent=2))
if __name__=='__main__':main()
