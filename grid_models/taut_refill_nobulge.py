# same test restricted to galaxies without a bulge (bulge mass-to-light is the least certain input)
import numpy as np, glob
from scipy.optimize import minimize_scalar
exec(open("taut_refill.py").read().split("# ---------------- SPARC")[0])
kpc=3.0857e19; c,H0=2.998e8,67.5e3/3.0857e22; a0g=c*H0/6
D=[np.loadtxt(f,comments="#",ndmin=2) for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat"))]
D=[d for d in D if np.all(d[:,5]==0)]
for Yd in (0.5,0.6):
    GO,GB,EL=[],[],[]
    for d in D:
        r,V,eV,Vg,Vd,Vb,_,_=d.T; gb=(Vg*abs(Vg)+Yd*Vd**2)*1e6/(r*kpc); go=V**2*1e6/(r*kpc)
        ok=(V>0)&(eV/V<0.1)&(gb>0)&(r>0); GO+=list(go[ok]);GB+=list(gb[ok]);EL+=list(2*eV[ok]/V[ok]/np.log(10))
    go,gb,el=map(np.array,(GO,GB,EL))
    print(f"\nbulge-free: {len(D)} galaxies, {len(go)} pts, M/L {Yd}")
    for lab,nu in nus.items():
        chi=lambda a0: np.sum(((np.log10(go)-np.log10(gb*nu(gb/a0)))/np.hypot(el,0.1))**2)
        b=minimize_scalar(lambda x: chi(10**x),bounds=(-10.4,-9.4),method="bounded")
        print(f"   {lab:38s} a0={10**b.x:.3e} chi2={b.fun:7.1f}   at cH0/6: {chi(a0g):7.1f}")
