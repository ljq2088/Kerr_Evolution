"""Read completed fresh a=.88 reports; compare source/response without fitting."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from environment_radial import RadialGreen
from environment_response import SampledResponse
from report_li_fixed_source_background import OriginalSourceCoordinate
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'docs/root_cause_followup_20260917'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def enc(z):return [float(np.real(z)),float(np.imag(z))]
def main():
    phasepath=OUT/'fresh_a088_cloud_phase_convention.json';pd=json.loads(phasepath.read_text())
    phase=complex(*pd['new_to_old_constant_phase_factor']);inputs={str(phasepath.relative_to(ROOT)):sha(phasepath)}
    fig,axes=plt.subplots(2,2,figsize=(11.6,7.4),sharex='col')
    rows=[]
    for j,orbit in enumerate([41.1,42.1]):
        f=OUT/f'fresh_a088_rp{str(orbit).replace(".","p")}_sl2_sm2_L6_q12_nr8_h32.json'
        d=json.loads(f.read_text());bp=ROOT/d['baseline_control']['file'];b=json.loads(bp.read_text())
        if d['status']!='truncated_single_mode_not_converged':raise ValueError('Fresh run incomplete')
        inputs[str(f.relative_to(ROOT))]=sha(f);inputs[str(bp.relative_to(ROOT))]=sha(bp)
        nf=SampledResponse.from_report(d);ob=SampledResponse.from_report(b)
        p=d['parameters'];op=b['parameters'];wc=d['background_provenance']['cloud']['temporal_omega'][0]
        w=wc+1/(orbit**1.5+.88)
        g=RadialGreen(.88,.3,w,2,2,rmax=32000,offset=1e-4,rtol=1e-11,infinity_method='coulomb')
        oldrp=1+np.sqrt(1-op['metric']['a']**2)
        fixed=SampledResponse(OriginalSourceCoordinate(g,oldrp),op['source_panels'],
            [s['r'] for s in b['samples']],[complex(*s['source']) for s in b['samples']],log_first=True)
        r=np.linspace(20.,200.,721);old=ob.evaluate(r)[0];fresh=phase*nf.evaluate(r)[0];control=fixed.evaluate(r)[0]
        for v,label,style in [(old,'Synchronized cloud + metric','-'),(control,'Old source, exact a=0.88 Green','--'),(fresh,'Fresh a=0.88 cloud + metric',':')]:
            axes[0,j].plot(r,1000*v.real,style,label=label,lw=1.8)
        axes[0,j].axhline(0,color='.7',lw=.6);axes[0,j].set_title(f'Orbit {orbit:g} M')
        scale=np.max(abs(control));diff=fresh-control
        axes[1,j].plot(r,100*diff.real/scale,label='Real difference',lw=1.5)
        axes[1,j].plot(r,100*diff.imag/scale,label='Imaginary difference',lw=1.5)
        axes[1,j].axhline(0,color='.7',lw=.6);axes[1,j].set_xlabel('Boyer–Lindquist r / M')
        axes[1,j].set_title('Fresh source effect at the same exact Green')
        refs={s['r']:s for s in b['samples']};pairs=[(s,refs[s['r']]) for s in d['samples'] if s['r'] in refs]
        sw=np.array([s['weight'] for s,_ in pairs]);nJ=np.array([phase*complex(*s['source']) for s,_ in pairs]);oJ=np.array([complex(*s['source']) for _,s in pairs])
        source_l2=float(np.sqrt(np.sum(sw*abs(nJ-oJ)**2)/np.sum(sw*abs(oJ)**2)))
        probe=np.array([20.,50.,100.,150.,200.]);nval=phase*nf.evaluate(probe)[0];oval=ob.evaluate(probe)[0];cval=fixed.evaluate(probe)[0]
        cp=OUT/f'fresh_a088_rp{str(orbit).replace(".","p")}_scalar22_comparison.json';c=json.loads(cp.read_text())
        inputs[str(cp.relative_to(ROOT))]=sha(cp)
        rows.append(dict(orbit=orbit,fresh_parameters=p,cloud_spectral_omega=d['background_provenance']['cloud']['spectral_omega'],
            fresh_source_common_nodes=len(pairs),common_nodes_lower_radius=min(s['r'] for s,_ in pairs),
            phase_aligned_source_weighted_relative_L2=source_l2,
            excluded_source_comparison='First32 logarithmic nodes differ because rplus differs; only64 identical outer nodes enter this L2.',
            raw_z_h=d['z_h'],phase_aligned_z_h=enc(phase*complex(*d['z_h'])),baseline_z_h=b['z_h'],
            raw_z_inf=d['z_inf'],phase_aligned_z_inf=enc(phase*complex(*d['z_inf'])),baseline_z_inf=b['z_inf'],
            raw_up_coefficient=d['up_coefficient'],baseline_up_coefficient=b['up_coefficient'],
            up_coefficient_caution='Bound Up endpoint normalization changes with kappa; compare physical fields.',
            fresh_flux=d['flux'],baseline_flux=b['flux'],
            fresh_horizon_flux_ratio=d['flux']['horizon']['orbital_energy']/b['flux']['horizon']['orbital_energy'],
            fresh_infinity_flux_ratio=(d['flux']['infinity']['orbital_energy']/b['flux']['infinity']['orbital_energy'] if b['flux']['infinity']['orbital_energy'] else None),
            samples=[dict(r=float(x),fresh_phase_aligned=enc(n),baseline=enc(o),fixed_source_exact_Green=enc(k),
                fresh_over_fixed_source_minus1=enc(n/k-1),fresh_over_baseline_amplitude=abs(n/o)) for x,n,o,k in zip(probe,nval,oval,cval)],
            grid_relative_L2_fresh_minus_fixed=float(np.linalg.norm(diff)/np.linalg.norm(control)),
            grid_max_absolute_difference_over_fixed_peak=float(np.max(abs(diff))/scale),
            minima=c['fresh']['absolute_field_local_minima'],baseline_minima=c['baseline']['absolute_field_local_minima'],
            wronskian_relative_spread=d['wronskian_relative_spread'],metric_sampling=d['metric_sampling'],
            normalization_audit=d['cloud_normalization_audit']))
    axes[0,0].set_ylabel(r'$10^3\,\mathrm{Re}\,R_{22}$')
    axes[1,0].set_ylabel('Difference / control peak (%)')
    axes[0,0].legend(fontsize=8);axes[1,0].legend(fontsize=8)
    for ax in axes.ravel():ax.grid(alpha=.2)
    fig.suptitle('Fresh fixed-spin scalar22 control; cloud phase fixed only from R_cloud(20)')
    fig.tight_layout();fig.savefig(OUT/'fresh_a088_scalar22_control.png',dpi=180);fig.savefig(OUT/'fresh_a088_scalar22_control.pdf')
    result=dict(status='finite_resolution_fresh_scalar22_completed',rows=rows,input_sha256=inputs,
        phase_convention=pd,no_fitted_amplitude_or_phase=True,
        implementation_sha256={str(Path(__file__).relative_to(ROOT)):sha(__file__),
            'src/report_li_fixed_source_background.py':sha(ROOT/'src/report_li_fixed_source_background.py')},
        limitations=['Full18 field has not been recomputed at a=.88.',
         'Fresh source and all metric entries use a=.88; the fixed-J control remains a deliberately inconsistent diagnostic.',
         'No claim of complete angular/radial convergence at L6/q12/nr8/h32.',
         'Global source phase is fixed from the independently computed cloud R(20), never fitted to a field or paper figure.',
         'Frozen temporal frequency and complex spectral spatial profile remain an explicitly approximate cloud; see cloud residual report.'])
    (OUT/'fresh_a088_scalar22_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    text=['# 固定 a=0.88 的新云与新源：scalar22 受控重算','',
        '两个轨道均已从新度规和新复谱云重新计算96个源节点。未复用旧源、未改通量前因子，也未调整参数拟合论文。这里只重算 scalar22；其余17个模式尚不是这个新背景。','',
        '统一参数：alpha=0.3，a=0.88，云频率 Re(omega)=0.2962935347256115，Im(omega)=2.216609396445795e-9。空间轮廓保留复谱；KS T=0 归一化后显式冻结时间增长。度规L6、角求积12、径向nr8/h32，源区间[rplus+0.0005,320]，受迫Green取Coulomb32000、视界offset0.0001。设置由旧文件实际metadata读取。','',
        '云R(20)整体相位为 -0.011727461259428496 rad。图与下方复场比较只乘单一因子 exp(+0.011727461259428496 i)，使它与旧同步云R(20)>0的约定一致；原始结果文件保留未旋转复值。该转换不改变任何通量、模长或径向模长谷。','',
        '| 轨道 | 新 FH22 | 旧 FH22 | 新/旧 FH22 | 新 FI22 | 旧 FI22 |','|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        h=row['fresh_flux']['horizon']['orbital_energy'];oh=row['baseline_flux']['horizon']['orbital_energy'];i=row['fresh_flux']['infinity']['orbital_energy'];oi=row['baseline_flux']['infinity']['orbital_energy']
        text.append(f"| {row['orbit']} | {h:.10e} | {oh:.10e} | {h/oh:.8f} | {i:.10e} | {oi:.10e} |")
    text+=['','通量为每 q^2(Mc/M) 的轨道能量通量，未乘 alpha^-3。束缚通道 FI=0；其 Up 系数依赖任意有限外端单位归一化，不能用系数大小比较物理场。','',
        '| 轨道 | 新源/旧源差：共同64节点加权L2 | 新源响应/固定旧源响应差：r20..200 L2 |','|---:|---:|---:|']
    for row in rows:text.append(f"| {row['orbit']} | {row['phase_aligned_source_weighted_relative_L2']:.6e} | {row['grid_relative_L2_fresh_minus_fixed']:.6e} |")
    text+=['','上表源比较仅使用 r>3 的64个完全相同节点；首32个对数节点随rplus变化，不纳入同点比较。响应差比较则使用相同的精确 a=0.88、精确实部频率 Green 核，区分源变化与近阈值传播变化。','']
    bound=next(x for x in rows if x['orbit']==42.1)
    text+=['rp42.1 的前三个径向模长谷：','', '| 背景 | 第一谷 r/M | 第二谷 r/M | 第三谷 r/M |','|---|---:|---:|---:|']
    for label,key in [('旧同步背景','baseline_minima'),('完整新a=0.88源','minima')]:
        text.append('| '+label+' | '+' | '.join(f"{v['r']:.8f}" for v in bound[key][:3])+' |')
    text+=['','这些谷与固定旧源但更换精确Green的诊断基本重合，说明新云/新度规源并未额外移动那些径向节点。不能将近零点处很大的相对场差解释为整体振幅同等改变。','',
        '质量切片核对：单位谱KS质量为1.000000000028788；冻结场同一切片能量为1.0000000269882492。有限BL切片在inner5e-4、outer1488.98给出0.9999998871498534；inner1e-6给出0.9999998039389256。sourceouter320时为0.9999998854658257。这些差异被原样记录、未吸收进新振幅，不能解释百分数级场偏差；有限BL切片也不是未来视界正则切片。','',
        '复谱云推导、局部KG缺陷、通量平衡和10项测试见 GENERAL_KERR_QUASIBOUND_CLOUD.md 与 general_kerr_quasibound_cloud_validation.json。本文计算没有宣布L6/q12/nr8/h32完全收敛，没有宣布完整论文场图已经复现。','',
        '文件：fresh_a088_scalar22_summary.json 保存输入SHA、两轨道通量、复场、源差、径向谷和归一化信息；fresh_a088_scalar22_control.png/.pdf 为不拟合相位/幅度的可复查图；两个 fresh_a088_rp*_sl2_sm2_L6_q12_nr8_h32.json 是完整标准samples布局。']
    (OUT/'FRESH_A088_SCALAR22_CONTROL.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print(json.dumps([dict(orbit=r['orbit'],horizon_ratio=r['fresh_horizon_flux_ratio'],infinity_ratio=r['fresh_infinity_flux_ratio'],source_L2=r['phase_aligned_source_weighted_relative_L2'],response_L2=r['grid_relative_L2_fresh_minus_fixed']) for r in rows],indent=2))
if __name__=='__main__':main()
