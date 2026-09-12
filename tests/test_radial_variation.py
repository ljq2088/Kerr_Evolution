"""Independent finite-difference checks of auxiliary-mass tangent solutions."""
import numpy as np
import pytest
from environment_radial import RadialGreen
from environment_radial_variation import RadialMassVariation


@pytest.mark.parametrize("a,omega,ell,m", [(.6,.1,1,1),(.877,.03,3,1),(.6,-.1,2,-1)])
def test_mass_tangent_against_independent_radial_solver(a,omega,ell,m):
    tangent=RadialMassVariation(a,omega,ell,m,rmax=300)
    step=.001*omega**2
    references=[RadialGreen(a,0,omega,ell,m,rmax=300,offset=1e-4,
                rtol=1e-12,mass_squared=x*step) for x in (-1,-.5,0,.5,1)]
    radii=np.geomspace(tangent.rmin+.01,300,40)
    for name in ("insol","upsol"):
        values=[getattr(g,name).sol(radii) for g in references]
        derivative=(4*(values[3]-values[1])/step
                    -(values[4]-values[0])/(2*step))/3
        actual=getattr(tangent,name).sol(radii)
        assert np.max(abs(actual[:2]-values[2]))/np.max(abs(values[2]))<2e-10
        assert np.max(abs(actual[2:]-derivative))/np.max(abs(derivative))<2e-8
    wronskians=tangent.wronskian(radii)
    for row in wronskians:
        assert np.max(abs(row-row[0]))/abs(row[0])<1e-9


@pytest.mark.parametrize("kwargs", [{"omega":0.}, {"mass_squared":.02}, {"offset":0.}])
def test_invalid_variation_domain(kwargs):
    parameters=dict(a=.6,omega=.1,ell=1,m=1,rmax=300)
    parameters.update(kwargs)
    with pytest.raises(ValueError):
        RadialMassVariation(**parameters)
