import pytest
from li_alignment_contract import alignment_issues,assert_comparison_ready,li_units,frequencies,REQUIRED

def test_old_threshold_and_grid_source_cannot_be_relabelled():
    run=dict(REQUIRED,a_over_M=.8771530275949366,source_method='2d_quadrature',mass_normalization='KS_T0')
    fields={x['field'] for x in alignment_issues(run)}
    assert {'a_over_M','source_method','mass_normalization','frequency_treatment'}<=fields
    with pytest.raises(ValueError):assert_comparison_ready(run)

def test_unknown_boundary_and_frequency_are_not_silently_accepted():
    run=dict(REQUIRED,cloud_state=[1,1,0])
    assert len(alignment_issues(run))==3

def test_declared_configuration_gate_does_not_prove_scientific_accuracy():
    run=dict(REQUIRED,cloud_state=[1,1,0])
    for k in ('frequency_treatment','mass_horizon_prescription','boundary_cutoffs'):
        run[k]={'value':'test fixture','evidence':'test fixture, not production evidence'}
    assert_comparison_ready(run)

def test_field_and_flux_units_are_distinct_powers():
    field,flux=li_units(1+2j,-3.,.3)
    assert field==pytest.approx((1+2j)/.027)
    assert flux==pytest.approx(-3/.000729)

def test_quadrupolar_frequency_uses_cloud_m_not_hardcoded_one():
    p=frequencies(.298451,2,0,20,.88)
    assert p['m_g']==-2
    assert p['omega']==pytest.approx(.298451-2*p['Omega_g'])
