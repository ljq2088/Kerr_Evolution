import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from lorenz_spin1 import spin1_amplitudes
from lorenz_chi import chi_amplitudes
from lorenz_ghp import KerrGHP
from lorenz_spin1 import cky_tensor


def test_cky_defining_equation():
    g=KerrGHP(6.,1.2,.6,order=4)
    f=cky_tensor(g)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                derivative=f[a][b].derivative(c-1) if c in (1,2) else f[a][b]*0
                derivative-=sum(g.gamma[k][c][a]*f[k][b]+g.gamma[k][c][b]*f[a][k] for k in range(4))
                target=g.g[b][c]*g.g[a][0]-g.g[a][c]*g.g[b][0]
                np.testing.assert_allclose(derivative.value,target.value,atol=1e-11)


def test_spin1_source_against_table_with_explicit_convention():
    omega=2/(6**1.5+.6)
    refs={-1:[-3.2451604624131902-6.232200824947201j,.4512743289756237+.1071192157797318j],
          1:[540.2135753511935+1037.4585568097918j,-.4252591398551278+.290016360574153j]}
    for spin in (-1,1):
        mapped=np.asarray(spin1_amplitudes(6.,spin=spin))*np.sqrt(2)/omega**2
        np.testing.assert_allclose(mapped,refs[spin],rtol=1e-10)


def test_compact_chi_source_against_table_with_explicit_convention():
    omega=2/(6**1.5+.6)
    mapped=np.asarray(chi_amplitudes(6.))*(-1j)/(2*omega)
    refs=[38.78733745778548-21.44877783639125j,-.3487514080653916-.44902453507200246j]
    np.testing.assert_allclose(mapped,refs,rtol=1e-10)
