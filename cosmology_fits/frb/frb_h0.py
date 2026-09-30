"""H0 from free electrons: the Macquart relation with 94 localized FRBs (FRBs/FRB repository).
DM_obs = DM_ISM (NE2001, from repo) + DM_MWhalo (50) + DM_cosmic + DM_host/(1+z)
<DM_cosmic>(z) = (3 c H0 Ω_b f_d χ_e)/(8π G m_p) ∫ (1+z)/E(z) dz,  Ω_b h^2 = 0.02242 fixed (BBN/CMB)  ->  ∝ f_d / h
P(Δ = DM_cosmic/<DM_cosmic>) : Macquart+2020 form (α=β=3, σ = F z^-1/2, F = 0.32);  DM_host lognormal (μ, σ_h) fitted.
Main systematic: f_d, the fraction of baryons in diffuse ionised gas (0.84 standard; 0.75-0.93 range)."""
import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
rows=np.load("frb_rows.npy",allow_pickle=True)
z=np.array([r[1] for r in rows],float); DM=np.array([r[3] for r in rows],float); ism=np.array([r[4] for r in rows],float)
ok=(z>0.01); z,DM,ism=z[ok],DM[ok],ism[ok]
import os
HALO=float(os.environ.get('HALO',50)); DMex=DM-ism-HALO
keep=DMex>20.; print('dropped (DM_ex<20):',(~keep).sum()); z,DM,ism,DMex=z[keep],DM[keep],ism[keep],DMex[keep]
print(f"{ok.sum()} FRBs (z>0.01), z median {np.median(z):.2f}")
c=2.998e8; G=6.674e-11; mp=1.6726e-27; Mpc=3.0857e22; pc_cm3=3.0857e16*1e6
wb=0.02242; chi=0.875
def Ez(zz,Om,beta):
    if beta==0: return np.sqrt(Om*(1+zz)**3+1-Om)
    return np.array([brentq(lambda e: e*e-Om*(1+x)**3-(1-Om)*(e/(1+x))**-beta,1e-3,1e6) for x in np.atleast_1d(zz)])
def meanDM(zz,h,fd,Om=0.31,beta=0):
    zg=np.linspace(0,zz.max(),400); E=Ez(zg,Om,beta); I=np.concatenate([[0],np.cumsum(0.5*((1+zg[1:])/E[1:]+(1+zg[:-1])/E[:-1])*np.diff(zg))])
    H100=100e3/Mpc; K=3*c*(H100**2)*wb/(8*np.pi*G*mp)/(h*H100)*fd*chi   # m^-2 per unit integral (c/H0 factor included)
    return np.interp(zz,zg,I)*K/pc_cm3*1e-6*1e6/1e6*1e6  # -> pc/cm^3
# unit check with standard numbers: <DM>(z=1) ~ 900-1000 pc/cm3 for h=0.7, fd=0.84
def Pcos_grid(sig):
    D=np.geomspace(0.01,60,2500); a=b=3.; lw=np.log(np.gradient(D))
    def lp(C): return -b*np.log(D)-(D**-a-C)**2/(2*a*a*sig*sig)+lw
    f=lambda C: np.exp(logsumexp(lp(C)+np.log(D))-logsumexp(lp(C)))-1
    C0=brentq(f,-30,30)
    l=lp(C0)-lw; p=np.exp(l-l.max()); p/=np.trapezoid(p,D); return D,p
F=float(os.environ.get('FEED',0.32)); PC=[Pcos_grid(F/np.sqrt(zi)) for zi in z]
mus=np.log(np.linspace(40,300,27)); shs=np.linspace(0.3,1.3,11)
Dg=PC[0][0]; lw=np.log(np.gradient(Dg)); LP=np.array([np.log(p+1e-300) for _,p in PC])+lw   # (N, nD)
def loglike_grid(h,fd,beta=0):
    m=meanDM(z,h,fd,beta=beta)
    host=(DMex[:,None]-m[:,None]*Dg[None,:])*(1+z[:,None])            # (N,nD)
    good=host>0.5; lx=np.log(np.where(good,host,1.))
    out=np.empty((len(mus),len(shs)))
    for i,mu in enumerate(mus):
        for j,sh in enumerate(shs):
            lh=np.where(good,-lx-np.log(sh*np.sqrt(2*np.pi))-(lx-mu)**2/(2*sh*sh)+np.log(1+z[:,None]),-np.inf)
            out[i,j]=np.sum(logsumexp(lh+LP,axis=1))
    return out
print(f"check <DM_cosmic>(z=1) = {meanDM(np.array([1.0]),0.7,0.84)[0]:.0f} pc/cm3 (literature ~ 950)")
hs=np.linspace(0.45,1.10,66)
def posterior(fd,beta=0):
    LL=np.array([loglike_grid(h,fd,beta) for h in hs])
    Lh=logsumexp(LL.reshape(len(hs),-1),axis=1); P=np.exp(Lh-Lh.max()); P/=P.sum()
    mean=np.sum(hs*P); cdf=np.cumsum(P); lo,med,hi=[np.interp(q,cdf,hs) for q in (0.16,0.5,0.84)]
    return med,lo,hi,P
import sys
mode=sys.argv[1]
if mode=="fd":
    for fd in [float(x) for x in sys.argv[2].split(",")]:
        med,lo,hi,P=posterior(fd)
        print(f"HALO={HALO} F={F} f_d={fd}: H0 = {100*med:.1f} (+{100*(hi-med):.1f} / -{100*(med-lo):.1f})  P(<69) {P[hs<0.69].sum():.2f}  P(>73) {P[hs>0.73].sum():.2f}", flush=True)
else:
    med,lo,hi,P=posterior(0.84,beta=0.5)
    print(f"β=½ law, HALO={HALO} F={F} f_d=0.84: H0 = {100*med:.1f} (+{100*(hi-med):.1f} / -{100*(med-lo):.1f})", flush=True)
