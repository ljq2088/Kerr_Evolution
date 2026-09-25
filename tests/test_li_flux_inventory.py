from report_li_aligned_flux import inventory,summarize

def test_dipole_flux_inventory():
    modes=inventory(1)
    assert len(modes)==18 and len(set(modes))==18
    assert (0,0) in modes and (2,2) in modes and (5,-5) in modes
    assert all(m!=1 and (ell+m)%2==0 for ell,m in modes)
    assert len([(ell,m) for ell,m in modes if m>=2])==6

def test_missing_mode_is_never_published_as_total():
    expected=inventory(1)
    partial=[dict(ell=0,m=0,flux={b:dict(orbital_energy=1.) for b in ('infinity','horizon')})]
    out=summarize(partial,expected,20.,1)
    assert not out['complete'] and 'total_flux' not in out and len(out['missing_modes'])==17
