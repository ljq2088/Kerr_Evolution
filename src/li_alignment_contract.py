"""Reject method/parameter mixing before claiming a Li-aligned comparison."""
import math

REQUIRED={
    'reference':'2507.02045v2', 'a_over_M':.88,'alpha':.3,
    'signature':'-+++','fourier':'exp(-i omega t+i m phi)',
    'angular_norm':'unit_sphere','cloud_radial_method':'Leaver',
    'cloud_radial_terms':150,'metric_basis':'weighted_tetrad_spin_spherical',
    'metric_spherical_lmax':18,'source_method':'semi_analytical_separated',
    'fourier_pmax':12,'fourier_dps':64,
    'horizon_series_order':4,'infinity_series_order':4,
    'mass_normalization':'BL_Killing_Li_zeta_def',
}

def alignment_issues(run):
    issues=[dict(field=k,expected=v,actual=run.get(k)) for k,v in REQUIRED.items() if run.get(k)!=v]
    for k in ('frequency_treatment','mass_horizon_prescription','boundary_cutoffs'):
        evidence=run.get(k)
        if not isinstance(evidence,dict) or not evidence.get('value') or not evidence.get('evidence'):
            issues.append(dict(field=k,expected='documented value and evidence',actual=evidence))
    if run.get('cloud_state') not in ([1,1,0],[2,2,0]):
        issues.append(dict(field='cloud_state',expected='[ell_c,m_c,n_c]=[1,1,0] or [2,2,0]',actual=run.get('cloud_state')))
    return issues

def assert_comparison_ready(run):
    issues=alignment_issues(run)
    if issues:raise ValueError('Li alignment incomplete: '+', '.join(x['field'] for x in issues))
    # This validates declared configuration, NOT underlying implementation or numerical accuracy.

def li_units(unit_mass_field,unit_mass_flux,alpha):
    if not math.isfinite(alpha) or alpha<=0:raise ValueError('alpha must be positive and finite')
    return unit_mass_field/alpha**3,unit_mass_flux/alpha**6

def frequencies(omega_c,m_c,m,r0,a):
    if r0<=1+math.sqrt(1-a*a):raise ValueError('Orbit outside horizon required')
    orbital=1/(r0**1.5+a);mg=m-m_c
    return dict(Omega_g=orbital,m_g=mg,omega_g=mg*orbital,omega=omega_c+mg*orbital)
