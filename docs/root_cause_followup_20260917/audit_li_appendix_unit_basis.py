"""Independent pointwise audit of Li v2 Appendix source formula.
Diagnostic only: no production physics modifications. Ten weighted unnormalized
Kinnersley components are unit-input fixtures, not physical point-particle modes.
Both incompatible printed L/Ldag mappings are preserved explicitly.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from environment_source import kerr_metric,connection,ThresholdCloud,angular_mode
OUT=Path(__file__).resolve().parent
NAMES=['lp_lp','lm_lm','mp_mp','mm_mm','Gamma_lp_mp','barGamma_lp_mm','barGamma_lm_mp','Gamma_lm_mm','SigmaDelta_lp_lm','trace']
def pair(z): return [float(np.real(z)),float(np.imag(z))]
def fixture(r,t,a,w,m):
    # Entire smooth separated field, without imposing KG: no on-shell cancellation
    # can hide an error in the Hessian contraction.
    beta=-.13+.021j
    R=np.exp(beta*r)*(1+.02*r*r)
    Rp=np.exp(beta*r)*(beta*(1+.02*r*r)+.04*r)
    Rpp=np.exp(beta*r)*(beta*beta*(1+.02*r*r)+.08*beta*r+.04)
    S=np.sin(t)+.17*np.cos(t)**2
    Sp=np.cos(t)-.34*np.sin(t)*np.cos(t)
    Spp=-np.sin(t)-.34*(np.cos(t)**2-np.sin(t)**2)
    f=R*S;grad=np.array([-1j*w*f,Rp*S,R*Sp,1j*m*f])
    partial=np.empty((4,4),complex)
    partial[0]=-1j*w*grad;partial[:,0]=partial[0]
    partial[3]=1j*m*grad;partial[:,3]=partial[3]
    partial[1,1]=Rpp*S;partial[2,2]=R*Spp
    partial[1,2]=partial[2,1]=Rp*Sp
    inv,christoffel=connection(r,t,a)
    H=partial-np.einsum('kij,k->ij',christoffel,grad)
    return R,Rp,Rpp,S,Sp,Spp,H,inv

def geometry(r,t,a):
    d=r*r-2*r+a*a;G=r+1j*a*np.cos(t);B=G.conjugate();sig=G*B
    U=np.array([[(r*r+a*a)/d,1,0,a/d],[-(r*r+a*a)/d,1,0,-a/d],
        [1j*a*np.sin(t),0,1,1j/np.sin(t)],[-1j*a*np.sin(t),0,1,-1j/np.sin(t)]])
    return d,G,B,sig,U

def covariant_from_weighted(x,r,t,a):
    d,G,B,sig,U=geometry(r,t,a);q=np.zeros((4,4),complex)
    for n in range(4):q[n,n]=x[n]
    for i,j,v in [(0,2,x[4]/G),(0,3,x[5]/B),(1,2,x[6]/B),(1,3,x[7]/G),(0,1,x[8]/(sig*d))]:
        q[i,j]=q[j,i]=v
    q[2,3]=q[3,2]=sig*x[9]-d*q[0,1]
    ui=np.linalg.inv(U)
    return ui@q@ui.T,q,U

def appendix(x,r,t,a,w,m,F,angular_mapping):
    R,Rp,Rpp,S,Sp,Spp=F[:6]
    d,G,B,sig,_=geometry(r,t,a);st=np.sin(t);ct=np.cos(t);cot=ct/st
    K=(r*r+a*a)*w-a*m;V=K/d;Vp=(2*r*w*d-K*(2*r-2))/d**2
    D=Rp-1j*V*R;Dd=Rp+1j*V*R
    DD=Rpp-2j*V*Rp+(-1j*Vp-V*V)*R
    DdDd=Rpp+2j*V*Rp+(1j*Vp-V*V)*R
    # Appendix map delta=Ldag/(sqrt2 Gamma) fixes L to m_minus.
    # Body line 423 instead explicitly assigns L to m_plus.
    sign=1 if angular_mapping=='appendix_delta' else -1
    Q=sign*(m/st-a*w*st);Qp=sign*(-m*ct/st**2-a*w*ct)
    L=Sp+Q*S;Ld=Sp-Q*S
    Lm1L=Spp+(2*Q-cot)*Sp+(Qp+Q*Q-cot*Q)*S
    Ldm1Ld=Spp+(-2*Q-cot)*Sp+(-Qp+Q*Q+cot*Q)*S
    L1Ld_plus_Ld1L=2*(Spp+cot*Sp-Q*Q*S)
    s2=(x[2]*Lm1L+x[3]*Ldm1Ld)*R/(4*sig)
    s1=d/(2*G**3*B)*(x[4]*(1j*a*st*Dd*S-R*L)+x[7]*(1j*a*st*D*S-R*Ld))
    s1+=d/(2*G**2*B)*(x[4]*L*Dd+x[7]*Ld*D)
    s1+=d/(2*G*B**2)*(x[5]*Ld*Dd+x[6]*L*D)
    s1-=d/(2*G*B**3)*(x[5]*(1j*a*st*Dd*S+R*Ld)+x[6]*(1j*a*st*D*S+R*L))
    common=sig*sig*x[9]-2*x[8]
    s0=common/(4*G**3*B**2)*(d*(D+Dd)*S+1j*a*st*R*(L+Ld))
    s0+=common/(4*G**2*B**3)*(d*(D+Dd)*S-1j*a*st*R*(L+Ld))
    s0+=(sig*sig*x[9]-x[8])/(4*sig*sig)*R*L1Ld_plus_Ld1L
    s0+=d*d/(4*sig)*(x[1]*DD+x[0]*DdDd)*S
    return s2+s1+s0


def source_targets(x,r,t,a,F):
    H,inv=F[6:];sig=geometry(r,t,a)[3]
    h,q,U=covariant_from_weighted(x,r,t,a)
    direct=sig*np.einsum('ij,ij->',inv@h@inv,H)
    eta=U@kerr_metric(r,t,a)@U.T;eta_inv=np.linalg.inv(eta)
    tetrad=sig*np.einsum('ij,ij->',eta_inv@q@eta_inv,U@H@U.T)
    return direct,tetrad

def cloud_fixture(cloud,r,t):
    R,Rp=cloud.radial(r)[:,0];Rpp=cloud.rhs(r,[R,Rp])[1]
    S,Sp,A=angular_mode(t,1,1,cloud.c2)
    Spp=-np.cos(t)/np.sin(t)*Sp-(cloud.c2*np.cos(t)**2-1/np.sin(t)**2+A)*S
    _,H,inv=cloud.hessian(r,t)
    return R,Rp,Rpp,S,Sp,Spp,H,inv

def static_massless_fixture(r,t):
    # Exact Schwarzschild static ell=1,m=0 massless KG solution.
    # Delta R''+Delta' R' - 2R = 2(r-1)-2(r-1)=0.
    R,Rp,Rpp=r-1.,1.,0.;S,Sp,Spp=np.cos(t),-np.sin(t),-np.cos(t)
    grad=np.array([0.,Rp*S,R*Sp,0.],complex)
    partial=np.zeros((4,4),complex);partial[1,1]=Rpp*S;partial[2,2]=R*Spp
    partial[1,2]=partial[2,1]=Rp*Sp
    inv,gamma=connection(r,t,0.)
    return R,Rp,Rpp,S,Sp,Spp,partial-np.einsum('kij,k->ij',gamma,grad),inv

def omitted_radial_term(x,r,t,a,w,m,F):
    R,Rp,Rpp,S=F[:4];d,_,_,sig,_=geometry(r,t,a)
    K=(r*r+a*a)*w-a*m
    return x[8]*((d*Rpp+(2*r-2)*Rp+K*K/d*R)*S)/(2*sig*sig)

def main():
    fixtures=[]
    for case in [dict(a=0.,r=5.,theta=.71,w=.2963,m=1),dict(a=.6,r=3.1,theta=1.1,w=.23,m=1),
              dict(a=.8771530275949366,r=10.,theta=.8,w=.296287,m=1),
              dict(a=.88,r=45.,theta=2.1,w=.2962935347256115,m=1)]:
        fixtures.append((dict(kind='manufactured_off_shell',**case),fixture(case['r'],case['theta'],case['a'],case['w'],case['m'])))
    cloud=ThresholdCloud()
    for r,t in [(3.1,.71),(10.,.8),(45.,2.1)]:
        case=dict(kind='exact_synchronized_massive_cloud',a=cloud.a,r=r,theta=t,w=cloud.omega,m=1,mu=cloud.mu,
            angular_A=cloud.lam,lambda_c=float(cloud.lam+cloud.a**2*cloud.omega**2-2*cloud.a*cloud.omega))
        fixtures.append((case,cloud_fixture(cloud,r,t)))
    for r,t in [(3.,.71),(10.,1.1)]:
        fixtures.append((dict(kind='exact_static_massless_Schwarzschild_ell1',a=0.,r=r,theta=t,w=0.,m=0,mu=0.,angular_A=2.,lambda_c=2.),
                         static_massless_fixture(r,t)))
    rows=[];cloud_checks=[]
    for case,F in fixtures:
        a,r,t,w,m=[case[k] for k in ['a','r','theta','w','m']]
        H,inv=F[6:];sig=geometry(r,t,a)[3];field=F[0]*F[3]
        if 'mu' in case:
            box=np.einsum('ij,ij->',inv,H)
            cloud_checks.append(dict(**case,KG_residual=pair(box-case['mu']**2*field),
                field=pair(field),scalar_Hessian_trace=pair(box)))
        basis=np.eye(10,dtype=complex)
        point_scale=max(abs(source_targets(x,r,t,a,F)[0]) for x in basis)
        for j,name in enumerate(NAMES):
            x=basis[j]
            target,tetrad=source_targets(x,r,t,a,F)
            results={mapping:appendix(x,r,t,a,w,m,F,mapping) for mapping in ['appendix_delta','literal_body_423']}
            omitted=omitted_radial_term(x,r,t,a,w,m,F)
            correction_error=target-results['appendix_delta']-omitted
            row=dict(**case,component=name,coordinate_contraction=pair(target),tetrad_contraction=pair(tetrad),
                tetrad_absolute_error=float(abs(target-tetrad)),
                appendix={k:dict(value=pair(z),absolute_error=float(abs(z-target)),
                    relative_error=float(abs(z-target)/abs(target)) if abs(target)>1e-12*point_scale else None,
                    scaled_error=float(abs(z-target)/point_scale)) for k,z in results.items()},
                point_coefficient_norm=float(point_scale),
                omitted_radial_term=pair(omitted),analytic_difference_absolute_error=float(abs(correction_error)),
                analytic_difference_scaled_error=float(abs(correction_error)/point_scale))
            if 'mu' in case:
                on_shell=x[8]*(case['mu']**2*r*r+case['lambda_c'])*field/(2*sig*sig)
                row['on_shell_omitted_term']=pair(on_shell)
                row['on_shell_term_absolute_error']=float(abs(on_shell-omitted))
            rows.append(row)
    result=dict(status='pointwise_unit_basis_identity_audit_not_particle_source_or_flux_test',
        scalar_fixtures=['Entire manufactured separated complex function (off shell)',
            'Production exact synchronized ThresholdCloud |211> at three generic BL points',
            'Analytic exact massless Schwarzschild static ell1: (r-1)cos(theta), at two points'],
        target='Sigma*h^{mu nu}*Hessian_mu_nu',component_order=NAMES,
        source_definition='Li 2507.02045v2 main.tex source_components 861-912, h_decomp 378-396',
        mappings={'appendix_delta':'L along m_minus and Ldag along m_plus, implied by delta table lines 829-837 and explicit tetrad 355-369',
            'literal_body_423':'L along m_plus and Ldag along m_minus, literal statement in line 423'},
        identified_difference='direct - printed(Appendix mapping) = x8 * (Delta Rpp + Delta_prime Rp + Kc^2/Delta R) * S / (2 Sigma^2)',
        on_shell_difference='x8 * (mu^2 r^2 + lambda_c) * Phi_c / (2 Sigma^2), with x8=Sigma Delta h_lplus_lminus and lambda_c=A_c+a^2 omega_c^2-2am_c omega_c',
        causal_scope='Printed algebra only. Does not demonstrate what the author ran, nor identify actual figure mismatch. No source or field rescaling applied.',
        worst_relative_by_mapping={k:max(row['appendix'][k]['relative_error'] for row in rows if row['appendix'][k]['relative_error'] is not None) for k in ['appendix_delta','literal_body_423']},
        worst_printed_scaled_by_mapping={k:max(row['appendix'][k]['scaled_error'] for row in rows) for k in ['appendix_delta','literal_body_423']},
        relative_error_policy='relative_error=null for coefficient <=1e-12 of pointwise ten-coefficient infinity norm; scaled_error always uses that norm',
        worst_tetrad_absolute=max(row['tetrad_absolute_error'] for row in rows),
        worst_analytic_difference_scaled=max(row['analytic_difference_scaled_error'] for row in rows),
        worst_on_shell_term_absolute=max(row.get('on_shell_term_absolute_error',0.) for row in rows),
        on_shell_checks=cloud_checks,rows=rows)
    assert result['worst_analytic_difference_scaled'] < 1e-12
    for row in rows:
        if row['component']!='SigmaDelta_lp_lm':
            assert row['appendix']['appendix_delta']['scaled_error'] < 1e-12
    for row in cloud_checks:
        assert abs(complex(*row['KG_residual'])) < 1e-12*max(1.,abs(complex(*row['field'])))
    closed_forms=[]
    for row in rows:
        if row['kind']=='exact_static_massless_Schwarzschild_ell1' and row['component']=='SigmaDelta_lp_lm':
            r,t=row['r'],row['theta']
            correct=2*np.cos(t)/r**4;printed=(3-r)*np.cos(t)/r**4
            assert abs(complex(*row['coordinate_contraction'])-correct)<1e-13
            assert abs(complex(*row['appendix']['appendix_delta']['value'])-printed)<1e-13
            closed_forms.append(dict(r=r,theta=t,correct=correct,printed=printed))
    result['closed_form_Schwarzschild_test']=closed_forms
    result['regression_assertions']='passed: corrected ten-channel identity <=1e-12, other nine channels <=1e-12, exact KG checks, analytic Schwarzschild counterexample'
    files=['outputs/paper_original_reference/li_2507_02045v2/source/main.tex','src/environment_source.py','src/environment_cloud.py']
    result['source_sha256']={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}
    result['diagnostic_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'li_appendix_unit_basis_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['worst_relative_by_mapping','worst_tetrad_absolute','worst_analytic_difference_scaled','worst_on_shell_term_absolute']}))
    for name in NAMES:
        print(name,max(row['appendix']['appendix_delta']['scaled_error'] for row in rows if row['component']==name))
if __name__=='__main__':main()
