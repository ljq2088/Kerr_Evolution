"""Evidence-scaled H00 numerical budget; observed changes are not rigorous bounds."""
import hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from source_provenance import source_fingerprint,validate_saved_samples
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'docs/environment_reproduction'
INPUTS={}
def load(name):
    p=D/name;raw=p.read_bytes();INPUTS[name]=hashlib.sha256(raw).hexdigest();return json.loads(raw)
def z(pair):return complex(*pair)
def enc(v):return [v.real,v.imag]
def H(d):return d['flux']['horizon']['orbital_energy']
def comparison(a,b):
    za,zb=z(a['z_h']),z(b['z_h'])
    return dict(z_before=enc(za),z_after=enc(zb),delta_z=enc(zb-za),relative_complex_amplitude_change=abs(zb-za)/abs(za),flux_before=H(a),flux_after=H(b),relative_flux_change=H(b)/H(a)-1)
def main():
    coarse=load('remaining_horizon_L4_nr8_q10_h32_20260916.json');fine=load('remaining_horizon_L4_nr16_q10_h64_20260916.json')
    for d in [coarse,fine]:
        if 'z_h' not in d:raise ValueError('Fresh source calculation is incomplete')
        validate_saved_samples(d,source_fingerprint())
    assert coarse['source_provenance']==fine['source_provenance']
    changed={k:[v,fine['parameters'][k]] for k,v in coarse['parameters'].items() if v!=fine['parameters'][k]}
    assert set(changed)=={'radial_order','horizon_quadrature_order'}
    fresh=load('fresh_20260915_L18_mg-1_sl0.json')
    old4=load('forced_mode_nr8_nt10_L4_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    old4nr4=load('forced_mode_nr4_nt10_L4_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    old6=load('forced_mode_nr8_nt10_L6_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    old18=load('forced_mode_nr8_nt18_L18_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    l1q10=load('forced_mode_nr8_nt10_L1_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    l1q18=load('forced_mode_nr8_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h32.json')
    l1fine=load('forced_mode_nr16_nt18_L1_mg-1_sl0_inner0.0005_outer320_log_h64.json')
    coverage5=load('flux_coverage_L18_nt18_i6_h5_f5_m0go4000_m2go4000.json')
    coverage12=load('flux_coverage_L18_nt18_i12_h12_f12_m0go4000_m1go4000_m2go4000.json')
    assert coverage5['horizon']['computed_count']==coverage5['horizon']['required_count']
    assert coverage12['horizon']['computed_count']==coverage12['horizon']['required_count']
    assert all(coverage5['parameters'][k]==coverage12['parameters'][k] for k in ['alpha','rp','metric_ellmax','nr','nt','horizon_order','inner_offset','outer_source_cutoff'])
    total=coverage5['horizon']['finite_resolution_total'];total12=coverage12['horizon']['finite_resolution_total'];share=H(fresh)/total
    summed5={};summed12={}
    for dest,cover in [(summed5,coverage5),(summed12,coverage12)]:
        for m in cover['horizon']['modes']:dest[(m['ell'],m['m'])]=m['flux']['horizon']['orbital_energy']
    common_change=sum(summed12[key]-val for key,val in summed5.items())
    extra=sum(val for key,val in summed12.items() if key not in summed5)
    cutoff=load('horizon_source_cutoff_r20_audit.json')
    contact=load('orbit_source_matching_summary_20260916.json')
    independent=load('independent_scalar00_response_audit.json')
    reference=load('li_reference_comparison_20260916.json')
    ref=next(row for row in reference['rows'] if row['rp']==20 and row['boundary']=='horizon')
    old_guard=load('trace_kappa_domain_guard_validation.json')
    changed_sources={name:[digest,coarse['source_provenance']['files'].get(name)] for name,digest in fresh['source_provenance']['files'].items() if digest!=coarse['source_provenance']['files'].get(name)}
    historical=[dict(name='L4 old nr4 to old nr8; q10, h32',kind='historical_unversioned',**comparison(old4nr4,old4)),dict(name='L4 to L6; same nr8/q10/h32',kind='historical_unversioned',**comparison(old4,old6)),dict(name='L6/q10 to L18/q18; two controls changed',kind='historical_unversioned_combined',**comparison(old6,old18)),dict(name='L18 old to independently refreshed same grid',kind='historical_to_pre_guard_versioned',**comparison(old18,fresh)),dict(name='L1 q10 to q18; fixed radial grid',kind='historical_dipole_only',**comparison(l1q10,l1q18)),dict(name='L1 radial88 to176; q18',kind='historical_dipole_only',**comparison(l1q18,l1fine))]
    radial=comparison(coarse,fine);radial.update(changed_parameters=changed,sample_counts=[len(coarse['samples']),len(fine['samples'])],angular_order=10,metric_ellmax=4,source_fingerprint=coarse['source_provenance']['sha256'],source_hashes_match=True,current_provenance_validation='passed',limitation='Includes spin0/1 dipole and ell2..4 reconstruction. Not a full L18 new-grid convergence certificate.')
    projected_delta=z(fine['z_h'])-z(coarse['z_h']);projected_flux_change=abs((z(fresh['z_h'])+projected_delta)/z(fresh['z_h']))**2-1
    bridged=dict(old_unversioned_to_fresh_L4=comparison(old4,coarse),fresh_L4_to_pre_guard_fresh_L18=comparison(coarse,fresh),source_file_hash_changes_since_L18=changed_sources,guard_actual_point_equality=old_guard['trace_all_coefficients_array_equal'] and old_guard['kappa_all_coefficients_array_equal'],interpretation='This is an explicit historical bridge, not rebadging older source samples with current provenance. Guard equality was checked at one valid point; current fresh pair has its own identical provenance. Historical comparisons are not rigorous error bounds.',sensitivity_if_only_measured_L4_delta_is_transferred_to_full_L18=dict(delta_z=enc(projected_delta),relative_H00_flux_change=projected_flux_change,relative_total_flux_change=share*projected_flux_change))
    # Observed sensitivities, not independent stochastic errors or proven remainder bounds.
    endpoint=max(abs(r['contact_by_ell'][-1]['signed_flux_relative_change']) for r in contact['results'])
    outer=[]
    for row in fresh['outer_source_cutoff_sequence']:
        if row['outer_source_cutoff']>=80:
            outer.append(dict(cutoff=row['outer_source_cutoff'],z_h=row['z_h'],relative_flux_to_full=abs(z(row['z_h'])/z(fresh['z_h']))**2-1))
    e160=next(r for r in outer if r['cutoff']==160)['relative_flux_to_full']
    tail=max(abs(v) for v in cutoff['estimated_limit_relative_flux_change_range'])
    observations=[
        dict(label='New L4 radial 88 -> 176',scope='new controlled source pair; transfer only of measured L4 delta',relative_total_change=share*projected_flux_change),
        dict(label='Dipole inner source tail',scope='previous round; reused result, not recomputed; dipole only',relative_total_change=-share*tail),
        dict(label='L4 -> L6 metric cutoff',scope='historical same-grid H00 sensitivity',relative_total_change=share*historical[1]['relative_flux_change']),
        dict(label='L6/q10 -> L18/q18',scope='historical combined controls; not isolated angular/metric bound',relative_total_change=share*historical[2]['relative_flux_change']),
        dict(label='Scalar ell 5 -> 12 sum',scope='complete 18-channel vs85-channel historical sums; includes common-channel Green boundary changes',relative_total_change=(total12-total)/total),
        dict(label='Source outer 160 -> 320',scope='same pre-guard freshL18 source; 320->infinity not tested by this increment',relative_total_change=-share*e160),
        dict(label='Orbit contact projection',scope='maximum observed fullL18 contact diagnostic sensitivity; not a full bulk Lorenz norm bound',relative_total_change=share*endpoint),
        dict(label='Independent cloud/source/Green',scope='same metric, same source grid; shared-model errors excluded',relative_total_change=share*independent['horizon_flux_relative_change'])]
    modes=sorted([dict(ell=k[0],m=k[1],flux=v,fraction_total=v/total) for k,v in summed5.items()],key=lambda x:abs(x['flux']),reverse=True)
    required=[]
    for label,key in [('Updated Dyson','dyson_in_later_comparison_magnitude'),('Li independent','li_magnitude')]:
        target=-ref[key];needed=target-total;h00target=H(fresh)+needed
        required.append(dict(reference=label,target_total=target,local_relative_difference=total/target-1,needed_total_fraction=needed/total,needed_H00_flux_fraction_if_others_fixed=needed/H(fresh),needed_H00_amplitude_fraction_if_phase_unchanged=np.sqrt(h00target/H(fresh))-1))
    result=dict(status='controlled_low_ell_radial_test_and_evidence_budget_not_global_error_bound',parameters=dict(alpha=.3,a=coarse['parameters']['metric']['a'],rp=20,scalar_ell=0,scalar_m=0,metric_m=-1,inner_source_offset=.0005,outer_source_radius=320),new_radial_comparison=radial,historical_bridge=bridged,historical_convergence=historical,baseline=dict(total_horizon=total,fresh_H00=H(fresh),H00_fraction_total=share,coverage_identity='total uses historical complete coverage; fresh H00 differs only about6e-8 relative from its matching historical mode'),horizon_mode_weights=modes,scalar_cutoff=dict(ell5_total=total,ell12_total=total12,relative_total_change=(total12-total)/total,added_modes_flux=extra,common_modes_changed_flux=common_change,notes='Additional modes and changed Green outer radii are separated; scalar cutoff and metric reconstruction cutoff are distinct.'),source_outer_sequence=outer,observed_sensitivities=observations,observed_absolute_sum_not_error_bound=sum(abs(r['relative_total_change']) for r in observations),reference_comparison=required,reference_interpolation_sensitivities=dict(updated_dyson=ref['new_dyson_log_pchip_vs_linear_relative'],li=ref['li_log_pchip_vs_linear_relative'],interpretation='Difference between two interpolants on published vector curves, not a rigorous digitization uncertainty or numerical source error.'),unresolved=['No fresh fullL18 q18->q24 and radial88->176 combined source study exists at rp20.','No certified extrapolation of high-ell reconstruction cancellation into the inner source shell.','Published reference curves are not author numerical arrays; their radial sampling/interpolation uncertainty is not rigorously bounded.','Li uses a=.88 and a quasi-bound cloud frequency, while current calculation uses the exact synchronous a=.8771530275949366. Cloud imaginary-frequency treatment is not fully specified.','Independent module agreement cannot exclude shared source/matching/gauge/model assumptions.'],interpretation='Observed numerical sensitivities are much smaller than2-4%; their sum is only a scale comparison, not a confidence interval. Current evidence does not support attributing the remaining differences solely to demonstrated numerical error. Full paper reproduction remains incomplete.',inputs_sha256=INPUTS,implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out=D/'remaining_horizon_error_budget_20260916.json';out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    # Plot observed changes in percent, keeping the reference differences separate.
    fig,(ax,bx)=plt.subplots(1,2,figsize=(15,6),gridspec_kw={'width_ratios':[1,1.55]},constrained_layout=True)
    Ffine=H(fine);xs=[len(old4nr4['samples']),len(coarse['samples']),len(fine['samples'])];ys=[1e6*(H(t)/Ffine-1) for t in [old4nr4,coarse,fine]]
    ax.plot(xs,ys,'o-',color='#19558a',lw=1.7);ax.scatter(xs[0],ys[0],s=75,color='#9b9b9b',zorder=4,label='Historical nr4 point');ax.scatter(xs[1:],ys[1:],s=60,color='#19558a',zorder=5,label='Fresh controlled pair')
    ax.axhline(0,color='.65',lw=.7);ax.set(xlabel='Radial source nodes',ylabel='H00 flux difference from fresh fine value [ppm]',title='L4, q10: actual source quadrature');ax.legend(fontsize=9);ax.grid(alpha=.2)
    labels=[r['label'] for r in observations];vals=[100*abs(r['relative_total_change']) for r in observations]
    labels+=['Residual vs updated Dyson','Residual vs Li'];vals +=[100*abs(r['local_relative_difference']) for r in required]
    colors=['#19558a']+['#648aab']*(len(observations)-1)+['#bd4e35','#bd4e35']
    bx.barh(range(len(labels)),vals,color=colors);bx.set_yticks(range(len(labels)),labels,fontsize=9);bx.invert_yaxis();bx.set_xscale('log');bx.set_xlabel('Magnitude of relative total-horizon change [%]');bx.set_title('Observed sensitivities versus reference differences');bx.grid(axis='x',alpha=.25)
    for i,v in enumerate(vals):bx.text(v*1.18,i,f'{v:.3g}%',va='center',fontsize=8)
    bx.set_xlim(min(vals)*.4,max(vals)*6);fig.suptitle('Kerr alpha=0.3, rp=20: measured changes are not certified error bounds',fontsize=13)
    fig.savefig(D/'remaining_horizon_error_budget_20260916.png',dpi=180);fig.savefig(D/'remaining_horizon_error_budget_20260916.pdf');plt.close(fig)
    lines=['# Kerr 视界剩余偏差：数值误差预算复核','', '**目前不能认定剩余2%–4%只是数值误差。** 已实测的数值敏感性明显更小。本轮新增了包含非偶极项的实际源径向加密对照，未乘补偿因子、未改变物理公式。','', '## 本轮新增的受控对照','', '固定 alpha=.3、精确同步自旋、rp=20、scalar00、度规ell≤4、角向10节点、内外源边界和所有求解容差。只将各普通径向面板8点提高到16点、近视界对数面板32点提高到64点，即总源节点88→176。两次源均本轮真实计算，源码/包版本指纹一致。','',f'- 粗网格 Z_H = {z(coarse["z_h"]):.16g}；F_H = {H(coarse):.16g}。',f'- 精网格 Z_H = {z(fine["z_h"]):.16g}；F_H = {H(fine):.16g}。',f'- 复振幅相对变化 {radial["relative_complex_amplitude_change"]:.6g}；通量相对变化 {radial["relative_flux_change"]:.6g}（{100*radial["relative_flux_change"]:.6g}%）。', '', '这是包含spin0/1偶极以及ell2..4全部重构的直接源计算，补上此前只对L1做径向加密的缺口。它不是全L18新网格的误差认证。旧L4基准与本轮粗网格的变化单独记录；旧文件没有被伪装成当前版本。', '', '## 占比与已观察到的变化','', f'既有总视界通量中，scalar00占{100*share:.5f}%，scalar22约占{100*modes[1]["fraction_total"]:.5f}%。因此首先约束00通道，比只看极弱高阶模式的相对误差更有意义。','', '| 对照 | 相对总视界通量的变化幅值 | 证据范围 |','|---|---:|---|']
    for row in observations:lines.append(f'| {row["label"]} | {100*abs(row["relative_total_change"]):.6g}% | {row["scope"]} |')
    lines +=['',f'把上述变化幅值机械相加约为{100*result["observed_absolute_sum_not_error_bound"]:.6g}%。这个和只展示已测数值效应的量级：各项不是独立随机误差，且仍有未认证方向，**不能称为总误差上界或置信区间**。近视界内尾是复用上一轮已完成的结果，本轮没有重跑。','', f'标量求和从ell≤5扩至ell≤12的总变化为{100*(total12-total)/total:.6g}%；其中新增模式的通量是{extra:.8g}，共同模式因Green边界设置不同的变化是{common_change:.8g}。不能把标量输出阶数与度规重构阶数混为一谈。','', '## 剩余几个百分点意味着什么','']
    for row in required:lines.append(f'- 相对{row["reference"]}的差异为{100*row["local_relative_difference"]:+.4f}%。若其他模式不变，要靠00通道独自达到该参考值，其通量需改变{100*row["needed_H00_flux_fraction_if_others_fixed"]:+.4f}%，同相位幅度需改变{100*row["needed_H00_amplitude_fraction_if_phase_unchanged"]:+.4f}%。这些只是所需变化量，没有据此调整计算。')
    lines +=['','这些所需变化远大于本轮径向加密及已测截断敏感性。两条更新参考本身也不同，而且Li使用精确a=.88的准束缚云；本地使用a=.877153...的同步云。论文图形读数与参数/云处理差异不能作为本地离散误差直接相加。','', '## 仍缺少的认证','']
    lines +=['- '+s for s in result['unresolved']]
    lines +=['', '当前判断应是：大幅归一化偏差已解释；残余几百分点尚未归因。已测数值效应不足以支持“只是网格不够细”的结论，但尚不能排除未覆盖的重构/源构造误差或对比设置差异。下一步优先锁定同参数参考和云约定，再补完整L18角向与径向独立精化。','', '![数值敏感性与剩余偏差](environment_reproduction/remaining_horizon_error_budget_20260916.png)','', '机器结果保留每个输入的SHA256、源版本、完整复振幅及各项限制：`environment_reproduction/remaining_horizon_error_budget_20260916.json`。']
    (ROOT/'docs/REMAINING_HORIZON_ERROR_BUDGET_20260916.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(radial=radial,observed_sum_not_bound=result['observed_absolute_sum_not_error_bound'],reference=required),indent=2))
if __name__=='__main__':main()
