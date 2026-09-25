"""Generate separated monomials from the independently audited tetrad source.
Uses isolated SymPy in outputs/li_symbolic_dependencies; no production env change.
"""
import sys,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'outputs/li_symbolic_dependencies'))
import sympy as s
r,a,d,dp,V,Vp,w,m=s.symbols('r a d dp V Vp w m')
G,B,st,ct=s.symbols('Gamma barGamma sin_theta cos_theta')
R,Rp,Rpp,S,Sp,Spp=s.symbols('R Rp Rpp S Sp Spp')
x=s.symbols('x0:10');I=s.I;sig=G*B;cot=ct/st
D=Rp-I*V*R;Dd=Rp+I*V*R
DD=Rpp-2*I*V*Rp+(-I*Vp-V*V)*R
DdDd=Rpp+2*I*V*Rp+(I*Vp-V*V)*R
Q=m/st-a*w*st;Qp=-m*ct/st**2-a*w*ct
L=Sp+Q*S;Ld=Sp-Q*S
Lm1L=Spp+(2*Q-cot)*Sp+(Qp+Q*Q-cot*Q)*S
Ldm1Ld=Spp+(-2*Q-cot)*Sp+(-Qp+Q*Q+cot*Q)*S
mixed=2*(Spp+cot*Sp-Q*Q*S)
f=(x[2]*Lm1L+x[3]*Ldm1Ld)*R/(4*sig)
f+=d/(2*G**3*B)*(x[4]*(I*a*st*Dd*S-R*L)+x[7]*(I*a*st*D*S-R*Ld))
f+=d/(2*G**2*B)*(x[4]*L*Dd+x[7]*Ld*D)
f+=d/(2*G*B**2)*(x[5]*Ld*Dd+x[6]*L*D)
f-=d/(2*G*B**3)*(x[5]*(I*a*st*Dd*S+R*Ld)+x[6]*(I*a*st*D*S+R*L))
common=sig**2*x[9]-2*x[8]
f+=common/(4*G**3*B**2)*(d*(D+Dd)*S+I*a*st*R*(L+Ld))
f+=common/(4*G**2*B**3)*(d*(D+Dd)*S-I*a*st*R*(L+Ld))
f+=(sig**2*x[9]-x[8])/(4*sig**2)*R*mixed
f+=d*d/(4*sig)*(x[1]*DD+x[0]*DdDd)*S
correction=x[8]/(2*sig**2)*(d*Rpp+dp*Rp+d*V**2*R)*S
radial=(r,a,d,dp,V,Vp,w,m)
def encode(expr):
 rows=[]
 for term in s.Add.make_args(s.expand(expr)):
  powers=term.as_powers_dict()
  component=[i for i,v in enumerate(x) if powers.get(v,0)]
  ri=[i for i,v in enumerate((R,Rp,Rpp)) if powers.get(v,0)]
  si=[i for i,v in enumerate((S,Sp,Spp)) if powers.get(v,0)]
  assert len(component)==len(ri)==len(si)==1
  row=dict(component=component[0],radial_derivative=ri[0],angular_derivative=si[0],
      beta=int(-powers.get(G,0)),sigma=int(-powers.get(B,0)),
      sin_power=int(powers.get(st,0)),cos_power=int(powers.get(ct,0)),
      radial_powers=[int(powers.get(v,0)) for v in radial])
  coeff=term
  for var in (*radial,G,B,st,ct,R,Rp,Rpp,S,Sp,Spp,*x):coeff/=var**powers.get(var,0)
  assert not coeff.free_symbols
  row['coefficient']=[str(s.re(coeff)),str(s.im(coeff))]
  assert 0<=row['beta']<=3 and 0<=row['sigma']<=3
  rows.append(row)
 return rows
out=dict(status='source_operator_terms_not_full_pipeline',radial_variables=[str(v) for v in radial],
 angular_mapping='Appendix delta; fixed by tetrad directional derivative',
 printed_terms=encode(f),required_correction_terms=encode(correction))
(ROOT/'src/li_source_terms.json').write_text(json.dumps(out,indent=2)+'\n')
print('Terms',len(out['printed_terms']),'+',len(out['required_correction_terms']))
