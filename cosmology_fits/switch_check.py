import camb, numpy as np
from camb.baseconfig import camblib
import ctypes
camblib.set_vconst.argtypes=[ctypes.c_double]*2; camblib.set_vswitch.argtypes=[ctypes.c_double]*2
def run(me,zs,dz=30):
    camblib.set_vconst(1.0,me); camblib.set_vswitch(zs,dz)
    p=camb.set_params(H0=67.5,ombh2=0.0224,omch2=0.12,tau=0.055,As=2.1e-9,ns=0.965,lmax=2500)
    p.Recomb.use_rosenbrock=False
    r=camb.get_results(p); d=r.get_derived_params()
    return r.get_cmb_power_spectra(p,CMB_unit='muK')['total'][:,0], d['zstar'], d['thetastar']
ref,z0,t0=run(1.0,-1)
con,zc,tc=run(1.0078,-1)
print("const 1.0078 zstar",zc, "theta",tc, "(LCDM",z0,t0,")")
for zs in [5000,1300,1100,1000,900,800,600,300,150,20]:
    c,zz,tt=run(1.0078,zs)
    frac=np.sum(((c-ref)[30:2500])*(con-ref)[30:2500])/np.sum(((con-ref)[30:2500])**2)
    print(f"switch z={zs:5d}: zstar {zz:8.2f} theta {tt:.5f}  fraction of full effect {frac:5.2f}")
