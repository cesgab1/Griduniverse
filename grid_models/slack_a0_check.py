# exponential-slack law: chi2 with a0 fixed to cH0/6, and best a0 vs stellar mass-to-light
import numpy as np, glob
from scipy.optimize import brentq, minimize_scalar
kpc=3.0857e19; c,H0=2.998e8,67.5e3/3.0857e22; a0g=c*H0/6
def sig(e,a0): s0=a0/2; return e-s0*(1-np.exp(-e/s0))
def gof(gN,a0): return np.array([brentq(lambda g: sig(g,a0)-y,1e-20,y+10*a0) for y in gN])
D=[np.loadtxt(f,comments="#",ndmin=2) for f in sorted(glob.glob("/home/claude/sparc/r1/*_rotmod.dat"))]
for Yd in (0.4,0.5,0.6):
    Yb=1.4*Yd; GO=[];GB=[];EL=[]
    for d in D:
        r,V,eV,Vg,Vd,Vb,_,_=d.T; gb=(Vg*abs(Vg)+Yd*Vd**2+Yb*Vb**2)*1e6/(r*kpc); go=V**2*1e6/(r*kpc)
        ok=(V>0)&(eV/V<0.1)&(gb>0)&(r>0); GO+=list(go[ok]);GB+=list(gb[ok]);EL+=list(2*eV[ok]/V[ok]/np.log(10))
    go,gb,el=map(np.array,(GO,GB,EL))
    chi=lambda a0: np.sum(((np.log10(go)-np.log10(gof(gb,a0)))/np.hypot(el,0.1))**2)
    b=minimize_scalar(lambda x: chi(10**x),bounds=(-10.3,-9.6),method="bounded")
    print(f"Yd={Yd}: best a0={10**b.x:.3e}  chi2={b.fun:.1f} | a0=cH0/6={a0g:.3e}: chi2={chi(a0g):.1f}  (delta {chi(a0g)-b.fun:+.1f})")
