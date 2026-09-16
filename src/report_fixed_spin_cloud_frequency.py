"""Fixed-spin complex |211> Leaver spectrum; independent diagnostic only."""
import hashlib,json,time
from pathlib import Path
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
DPS=65
def pair(z,n=48):return [mp.nstr(mp.re(z),n),mp.nstr(mp.im(z),n)]
def cf(a,mu,w,terms,angular_terms=10):
    a,mu=mp.mpf(a),mp.mpf(mu);w=mp.mpc(w)
    rp=1+mp.sqrt(1-a*a);rm=2-rp;d=rp-rm
    q=-mp.sqrt(mu*mu-w*w)
    if mp.re(q)>=0:raise ValueError("Not the decaying infinity branch")
    beta=(mu*mu-2*w*w)/q-1;s=-1j*(2*rp*w-a)/d
    c2=a*a*(w*w-mu*mu)
    C=lambda l:mp.sqrt(mp.mpf(l*l-1)/(4*l*l-1)) if l>1 else mp.mpf(0)
    matrix=mp.matrix(angular_terms)
    for i in range(angular_terms):
        ell=1+2*i
        matrix[i,i]=ell*(ell+1)-c2*(C(ell+1)**2+C(ell)**2)
        if i<angular_terms-1:
            matrix[i,i+1]=matrix[i+1,i]=-c2*C(ell+1)*C(ell+2)
    values=mp.eig(matrix,left=False,right=False)
    A=min(values,key=lambda x:abs(x-2))
    sep=A+a*a*w*w-2*a*w
    def reduced(x):
        r=(rp-rm*x)/(1-x);delta=(r-rp)*(r-rm);K=(r*r+a*a)*w-a
        lp=s/x+beta/(1-x)+q*d/(1-x)**2
        lpp=-s/x**2+beta/(1-x)**2+2*q*d/(1-x)**3
        return (1-x)**2*(x*(lpp+lp*lp)+lp)+K*K/delta-mu*mu*r*r-sep
    x1,x2=mp.mpf(".25"),mp.mpf(".75")
    c1=(reduced(x2)-reduced(x1))/(x2-x1);c0=reduced(x1)-c1*x1
    defect=reduced(mp.mpf(".5"))-c0-c1/2
    if abs(defect)>mp.mpf("1e-50"):raise ValueError("Leaver reduction failed linearity")
    b0=1+2*s;b1=-2*b0+2*(q*d+beta);b2=b0-2*beta
    al=lambda n:(n+1)*(n+b0)
    be=lambda n:-2*n*(n-1)+b1*n+c0
    ga=lambda n:(n-1)*(n-2)+b2*(n-1)+c1
    ratio=mp.mpc(0)
    for n in range(terms,0,-1):
        ratio=-ga(n)/(be(n)+al(n)*ratio)
    return be(0)+al(0)*ratio,A,defect

def derived(a,wc,rp):
    Om=1/(rp**mp.mpf("1.5")+a)
    w=mp.re(wc)+Om
    kappa=mp.sqrt(mp.mpf(".3")**2-w*w)
    beta=2*w*w-mp.mpf(".3")**2
    return dict(a=mp.nstr(a,45),cloud_omega=pair(wc),
        orbital_omega=mp.nstr(Om,45),scalar_m=2,metric_m=1,
        scalar_real_omega=mp.nstr(w,45),mu_minus_scalar_real_omega=mp.nstr(mp.mpf(".3")-w,45),
        decay_kappa=mp.nstr(kappa,45),coulomb_beta=mp.nstr(beta,45),
        coulomb_eta=mp.nstr(beta/kappa,45),
        complex_decay_kappa=pair(mp.sqrt(mp.mpf(".3")**2-(wc+Om)**2)),
        frequency_prescription="Real cloud frequency for forced stationary response; complex kappa reported separately.")

def main():
    from paper_leaver_cloud import LeaverThresholdCloud
    output=ROOT/"docs/environment_reproduction/fixed_spin_cloud_frequency_20260917.json"
    with mp.workdps(DPS):
        a=mp.mpf(".88");mu=mp.mpf(".3")
        result=dict(status="fixed_spin_complex_Leaver_spectrum_running",
          physical_parameters=dict(a=".88",mu=".3",ell=1,m=1,radial_overtone=0),
          precision_decimal_digits=DPS,angular_terms=10,cases=[],
          recurrence="General x=(r-r+)/(r-r-) Frobenius reduction from paper_leaver_radial; minimal-solution backwards continued fraction from paper_leaver_cloud.",
          limitations=["Spectrum and derived threshold scales only; no cloud normalization, metric source, or forced radial field recomputed.",
            "Complex angular matrix is diagonalized rather than using the real symmetric eigensolver.",
            "Finite radial continued-fraction truncation is explicitly varied; numerical root residual alone is not truncation accuracy.",
            "No production module modified."])
        def save():
            temp=output.with_suffix(".tmp");temp.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n");temp.replace(output)
        save();seed=mp.mpc(".296294","1e-9")
        roots=[]
        for n in [150,300,500]:
            start=time.perf_counter()
            w=mp.findroot(lambda w:cf(a,mu,w,n)[0],(seed,seed+mp.mpf("1e-7")),tol=mp.mpf("1e-52"),maxsteps=40)
            residual,A,defect=cf(a,mu,w,n);roots.append(w);seed=w
            row=dict(radial_terms=n,omega=pair(w),angular_eigenvalue=pair(A),
                continued_fraction_residual=pair(residual),radial_reduction_linearity_residual=pair(defect),
                elapsed_seconds=time.perf_counter()-start)
            if len(roots)>1:row["change_from_previous_complex_omega"]=pair(w-roots[-2])
            result["cases"].append(row);save();print(json.dumps(row),flush=True)
        reference=LeaverThresholdCloud(mu=.3,terms=500,angular_terms=10,dps=DPS)
        residual,A,defect=cf(reference.a,mu,reference.omega,500)
        result["synchronous_generalization_check"]=dict(
            threshold_spin=mp.nstr(reference.a,48),threshold_omega=mp.nstr(reference.omega,48),
            specialized_cf_residual=mp.nstr(reference.residual,20),
            generalized_cf_residual=pair(residual),difference=pair(residual-reference.residual))
        w=roots[-1]
        result["fixed_spin_minus_stationary_cloud"]=dict(
            real_frequency_difference=mp.nstr(mp.re(w)-reference.omega,48),
            imaginary_frequency=mp.nstr(mp.im(w),48),
            spin_difference=mp.nstr(a-reference.a,48),
            fixed_spin_horizon_omega=mp.nstr(a/(2*(1+mp.sqrt(1-a*a))),48),
            rounding_error_using_0p296294=mp.nstr(mp.mpf(".296294")-mp.re(w),48))
        result["threshold_scales"]=[]
        for rp in [mp.mpf("42.1"),mp.mpf("41.8")]:
            baseline=derived(reference.a,reference.omega,rp)
            changed=derived(a,w,rp)
            onlyw=derived(reference.a,w,rp)
            rounded=derived(a,mp.mpf(".296294"),rp)
            row=dict(orbital_radius=mp.nstr(rp,10),stationary=baseline,fixed_spin_complex_spectrum=changed,
                change_cloud_frequency_only=onlyw,six_digit_cloud_frequency=rounded)
            for key in ["decay_kappa","coulomb_eta","scalar_real_omega"]:
                row[key+"_relative_change_fixed_spin_vs_stationary"]=mp.nstr(mp.mpf(changed[key])/mp.mpf(baseline[key])-1,40)
                row[key+"_relative_change_six_digit_vs_exact"]=mp.nstr(mp.mpf(rounded[key])/mp.mpf(changed[key])-1,40)
            result["threshold_scales"].append(row)
        result["implementation_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/"src/paper_leaver_cloud.py",ROOT/"src/paper_leaver_radial.py"]}
        result["status"]="fixed_spin_complex_Leaver_spectrum_completed"
        save();print(json.dumps(result["fixed_spin_minus_stationary_cloud"]),flush=True)
        print(json.dumps(result["threshold_scales"]),flush=True)
if __name__=="__main__":main()
