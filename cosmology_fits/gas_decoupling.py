# When does the cosmic gas stop following the CMB temperature, and does that set the electron switch?
import numpy as np, camb
from scipy.optimize import brentq
h=0.68; Om=0.31; OL=0.69; beta=0.63; Yp=0.245
p=camb.set_params(H0=100*h,ombh2=0.0224,omch2=Om*h*h-0.0224-0.00064,YHe=Yp)
r=camb.get_results(p)
z=np.logspace(0.5,3.3,4000)
ev=r.get_background_redshift_evolution(z,['x_e','T_b'],format='array')
xe,Tb=ev[:,0],ev[:,1]; Tc=2.7255*(1+z)
H=r.hubble_parameter(z)*1e3/3.0857e22            # 1/s
sT,me_c,arad=6.6524e-29,9.109e-31*2.998e8,7.5657e-16
fHe=Yp/(4*(1-Yp))
GC=8*sT*arad*Tc**4/(3*me_c)*xe/(1+fHe+xe)       # Compton heating rate of the gas
slope=np.gradient(np.log(Tb),np.log(1+z))     # 1 = tied to CMB, 2 = free adiabatic
def cross(y,target):
    i=np.where(np.diff(np.sign(y-target)))[0]; return [float(np.interp(target,[y[j],y[j+1]],[z[j],z[j+1]])) for j in i]
crit={'Compton rate = expansion rate':cross(GC/H,1.0),
      'gas 10% colder than CMB':cross(Tb/Tc,0.9),
      'halfway (T_gas ∝ (1+z)^1.5)':cross(slope,1.5),
      'gas 50% colder than CMB':cross(Tb/Tc,0.5)}
# energy budget: Δm/m needed if the switch is at z_s
rho_DE0=OL*1.0537e4*h*h; n_e0=0.0224*1.878e-29/1.6726e-24*(1-Yp/2)
def E(zz):
    a=1/(1+zz); e=np.sqrt(Om*a**-3+OL)
    for _ in range(50): e=np.sqrt(Om*a**-3+9.1e-5*a**-4+OL*(a*e)**-beta)
    return e
need_L=lambda zs: rho_DE0/(n_e0*(1+zs)**3*511e3)
need_B=lambda zs: rho_DE0*(E(zs)/(1+zs))**-beta/(n_e0*(1+zs)**3*511e3)
print("CMB fit: Δm/m = 0.78% ± 0.44%\n")
for k,v in crit.items():
    zs=max(v)
    print(f"{k:32s} z = {zs:6.0f}  -> predicted Δm/m: {need_L(zs):6.2%} (Λ), {need_B(zs):6.2%} (β law);  21-cm step at {1420.4/(1+zs):5.1f} MHz")
# per-electron (not gas) decoupling from photons
tC=3*me_c/(4*sT*arad*Tc**4); zi=cross(tC*H,1.0)
print(f"\nsingle free electron loses photon contact only at z ≈ {max(zi):.1f} -> needs Δm/m = {need_L(max(zi)):.0f}× (impossible)")
# coincidence odds: allowed window 5<z<900 (log-uniform); 1σ z_s range for Λ and β law
for lab,(lo,hi) in {'Λ':(134,205),'β law':(93,137)}.items():
    print(f"chance a random epoch lands in the {lab} 1σ range: {np.log((1+hi)/(1+lo))/np.log(901/6):.0%}")
np.savetxt('gas_decoupling.txt',np.c_[z,Tb/Tc,slope,GC/H],header='z Tgas/Tcmb slope Compton/H')

# Gradual version: the electron stays heavy in proportion to how tightly the gas is tied to the CMB, f = Γ_C/(Γ_C+H)
m=(z>12)                      # before reionisation
f=(GC/(GC+H))[m]; zz=z[m]
dF=np.gradient(f)             # f rises with z; released as f falls
w=(1+zz)**3*dF
for lab,fac in [('Λ',np.ones_like(zz)),('β law',np.array([(E(q)/(1+q))**beta for q in zz]))]:
    dm0=rho_DE0/(n_e0*511e3*np.sum(w*fac))
    zeff=(np.sum(w*fac)/np.sum(dF))**(1/3)-1
    print(f"gradual release ({lab}): predicted Δm/m = {dm0:.2%}, effective switch z ≈ {zeff:.0f}")
i=np.argsort(zz); print("release spans (10%-90%):", round(float(np.interp(0.1,f[i],zz[i]))), "-", round(float(np.interp(0.9,f[i],zz[i]))))
