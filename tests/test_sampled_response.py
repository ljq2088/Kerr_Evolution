import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from environment_response import SampledResponse
from environment_radial import RadialGreen
from environment_cloud import radial_coefficients


@pytest.mark.parametrize('omega',[.4,.2,-.4])
def test_continuous_manufactured_field_and_derivative(omega):
    g=RadialGreen(.7,.3,omega,2,2,rmax=100.)
    panels=[g.rmin,10.,20.,30.,80.]
    p=np.polynomial.Polynomial([0,0,0,0,1,-4,6,-4,1])*256
    def exact(r,derivative=0):
        x=(r-10)/20
        return np.where((x>0)&(x<1),p.deriv(derivative)(x)/20**derivative,0.)
    errors=[]
    for n in (8,16):
        radii=[]
        for i,(a,b) in enumerate(zip(panels[:-1],panels[1:])):
            x=np.polynomial.legendre.leggauss(n)[0]
            if i==0:
                a,b=np.log(a-g.rp),np.log(b-g.rp)
                radii.extend(g.rp+np.exp((a+b)/2+(b-a)*x/2))
            else:radii.extend((a+b)/2+(b-a)*x/2)
        r=np.array(radii)
        d,dp,v=radial_coefficients(r,.7,.3,omega,2,g.lam)
        J=d*exact(r,2)+dp*exact(r,1)+v*exact(r)
        response=SampledResponse(g,panels,r,J,log_first=True)
        probes=np.linspace(5,45,301)
        field,derivative=response.evaluate(probes)
        error=max(np.max(abs(field-exact(probes))),np.max(abs(derivative-exact(probes,1))))
        errors.append(error)
        # The source is smooth across the artificial panel boundary at 20.
        sides=response.evaluate([20-1e-7,20+1e-7])
        assert np.max(abs(sides[:,1]-sides[:,0]))<1e-6
    assert errors[-1]<2e-7
    assert errors[-1]<errors[0]/10


def test_incomplete_report_is_rejected():
    with pytest.raises(ValueError,match='completed'):
        SampledResponse.from_report(dict(status='source_sampling_in_progress'))


def test_known_carrier_recovers_rapidly_oscillating_manufactured_source():
    g=RadialGreen(.7,.3,.4,2,2,rmax=100.)
    panels=[10.,30.];x=np.polynomial.legendre.leggauss(20)[0]
    r=20+10*x;q=3.
    p=np.polynomial.Polynomial([0,0,0,0,1,-4,6,-4,1])*256
    z=(r-10)/20;carrier=np.exp(1j*q*r)
    f=p(z)*carrier
    df=(p.deriv()(z)/20+1j*q*p(z))*carrier
    ddf=(p.deriv(2)(z)/400+2j*q*p.deriv()(z)/20-q*q*p(z))*carrier
    d,dp,v=radial_coefficients(r,.7,.3,.4,2,g.lam)
    J=d*ddf+dp*df+v*f
    probe=np.linspace(10,30,301)
    exact=p((probe-10)/20)*np.exp(1j*q*probe)
    raw=SampledResponse(g,panels,r,J).evaluate(probe)[0]
    demodulated=SampledResponse(g,panels,r,J,phase=lambda radius:q*radius).evaluate(probe)[0]
    assert np.max(abs(raw-exact))>1e-3
    assert np.max(abs(demodulated-exact))<2e-6
