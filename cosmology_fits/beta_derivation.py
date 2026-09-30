# Riptide exchange: each grid link is crossed back and forth by returning flows; each crossing
# adds +1 or -1 to the link's tension imbalance. Crossings per stretch-time: N = c/(2 adot l) ∝ 1/adot.
# Net tension ∝ |sum| . If exchanges are independent -> ∝ sqrt(N) ∝ adot^-1/2 (β=1/2);
# fully coherent -> ∝ N (β=1). Test: exponent for exchanges with finite memory (correlation length Lc).
import numpy as np
rng=np.random.default_rng(1)
def net_rms(N,rho,trials=400):
    # AR(1) signs with correlation rho between successive exchanges
    s=np.empty((trials,N)); s[:,0]=rng.standard_normal(trials)
    e=rng.standard_normal((trials,N))*np.sqrt(1-rho**2)
    for i in range(1,N): s[:,i]=rho*s[:,i-1]+e[:,i]
    return np.sqrt(np.mean(np.sum(np.sign(s),axis=1)**2))
Ns=np.array([30,100,300,1000,3000,10000])
print("memory (corr. length)   local exponent d ln(tension)/d ln N  at N = 30..10000   -> β")
for rho in [0.0,0.5,0.9,0.99,0.999]:
    r=np.array([net_rms(n,rho,300 if n<3000 else 120) for n in Ns])
    ex=np.diff(np.log(r))/np.diff(np.log(Ns))
    Lc=-1/np.log(rho) if rho>0 else 0
    print(f"corr length {Lc:7.1f}    "+"  ".join(f"{x:.2f}" for x in ex))
# how large is N physically? comoving link = Planck length ... galaxy scale; adot today = H0
c=2.998e8; H0=2.2e-18
for lab,l in [('Planck length',1.6e-35),('proton',1e-15),('1 metre',1.0),('1 Mpc',3.1e22)]:
    print(f"link {lab:14s}: N ≈ {c/(2*H0*l):.1e} exchanges per stretch-time")
