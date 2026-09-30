# "Fast stretch" dark ingredient: pressureless on average (acts like dark matter in the background) but with a sound speed c_s,
# so it cannot settle in wells shallower than ~c_s. Linear growth from z=1000 (dark fluid + baryons), vs cold dark matter.
import numpy as np, camb
from scipy.integrate import solve_ivp
h=0.677; Om=0.311; Ob=0.049; Ox=Om-Ob; OL=1-Om
p=camb.set_params(H0=100*h,ombh2=Ob*h*h,omch2=Ox*h*h,As=2.1e-9,ns=0.965); p.set_matter_power(redshifts=[0],kmax=20)
r=camb.get_results(p); kh,_,P=r.get_matter_power_spectrum(minkh=1e-3,maxkh=10,npoints=400); P=P[0]
cH0=2997.9/h/h*h  # c/H0 in Mpc/h -> 2997.9 Mpc/h
def growth(k, cs):   # k in h/Mpc, cs in km/s
    q=(cs/2.998e5)**2*(k*2997.9)**2          # (c_s k / H0)^2
    def f(lna,y):
        a=np.exp(lna); E2=Om*a**-3+OL; dlnE=-1.5*Om*a**-3/E2
        dx,vx,db,vb=y; src=1.5*(Ox*a**-3*dx+Ob*a**-3*db)/E2
        return [vx, -(2+dlnE)*vx+src-q/(a*a*E2)*dx, vb, -(2+dlnE)*vb+src]
    a0=1e-3; y0=[a0,a0,0.2*a0,0.2*a0]      # baryons still catching up after decoupling
    s=solve_ivp(f,[np.log(a0),0],y0,rtol=1e-7,atol=1e-12)
    dx,db=s.y[0,-1],s.y[2,-1]; return (Ox*dx+Ob*db)/Om
def sig8(Pk):
    x=kh*8; W=3*(np.sin(x)-x*np.cos(x))/x**3; return np.sqrt(np.trapezoid(Pk*W**2*kh**2,kh)/(2*np.pi**2))
ref=np.array([growth(k,0.0) for k in kh])
s8=sig8(P)
print(f"ΛCDM σ8 = {s8:.3f}")
print("c_s (km/s)  c_s²/c²     σ8     P(k)/P_CDM at k=0.05  0.1  0.2  0.5 h/Mpc   Jeans scale today (Mpc/h)")
for cs in [50,100,150,240,300,540,1000]:
    g=np.array([growth(k,cs) for k in kh]); T2=(g/ref)**2
    kJ=np.sqrt(1.5*Ox)*100/cs        # h/Mpc
    print(f"{cs:6d}     {(cs/2.998e5)**2:8.1e}   {sig8(P*T2):.3f}     "+"  ".join(f"{np.interp(k,kh,T2):.2f}" for k in [0.05,0.1,0.2,0.5])+f"        {2*np.pi/kJ:6.0f}")
print("\nneeded window: galaxies (v_c 150-300 km/s) must NOT hold it -> c_s ≳ 300; clusters (σ ~ 1000) must -> c_s ≲ 1000")
print("limits on c_s²/c² (constant sound speed): CMB 3.2e-6 (99.7%; c_s < 540 km/s); CMB lensing + BAO + WiggleZ 0.2-0.63e-6 (95%; c_s < 130-240 km/s)")
print("measured σ8/S8 ≈ 0.77-0.83 (DES/KiDS/ACT lensing)")
