import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from report_nonstatic_paper_projection import transform, PAIRS
from environment_source import kerr_metric


def test_paper_tetrad_reconstructs_inverse_metric_and_background_projections():
    for r,t,a in ((4.,.7,0.),(6.,1.2,.6),(20.,2.4,.8771530275949366)):
        delta=r*r-2*r+a*a;sigma=r*r+a*a*np.cos(t)**2
        lp=np.array([(r*r+a*a)/delta,1.,0.,a/delta],complex)
        lm=np.array([-lp[0],1.,0.,-lp[3]])
        mp=np.array([1j*a*np.sin(t),0.,1.,1j/np.sin(t)]);mm=mp.conjugate()
        inverse=(delta*(np.outer(lp,lm)+np.outer(lm,lp))+np.outer(mp,mm)+np.outer(mm,mp))/(2*sigma)
        g=kerr_metric(r,t,a)
        np.testing.assert_allclose(inverse@ g,np.eye(4),atol=2e-14)
        weights=transform(r,t,a)
        h=np.array([g[i,j] for i,j in PAIRS])
        expected=np.zeros(10);expected[8]=2*sigma*sigma;expected[9]=4
        np.testing.assert_allclose(weights[0]@h,expected,atol=1e-9)
        step=1e-4
        numerical=(transform(r+step,t,a)[0]-transform(r-step,t,a)[0])/(2*step)
        np.testing.assert_allclose(weights[1],numerical,rtol=2e-8,atol=1e-8)
