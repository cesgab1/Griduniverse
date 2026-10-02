"""
Size of the predicted 21-cm 'step' (electron heavier by delta before the switch at z_s, relaxing to today's value after).
Standard dark-ages model (no stars): T_k from Compton coupling with residual ionisation 2e-4; H-H collisional coupling
kappa10 = 3.1e-11 T^0.357 exp(-32/T) cm^3/s; x_c = n_H kappa T*/(A10 T_g); T_s = (T_g + x_c T_k)/(1 + x_c);
T_b = 27 mK (1 - T_g/T_s) sqrt((1+z)/10) (wb/0.023) sqrt(0.15/wm).
Electron mass enters as: nu21 ~ m_e^2 (so T* ~ m_e^2), A10 ~ nu^3 mu_B^2 ~ m_e^4, optical depth ~ A10/nu^2 ~ unchanged.
Before the switch the line sits at nu21 (1 + 2 delta): observed frequency maps to a different redshift -> a jump at
nu_s = nu21/(1+z_s), plus the small change in coupling. Energy lost by electrons goes to the grid, not the gas (no heating).
"""
import numpy as np
from scipy.integrate import solve_ivp
wb, wm = 0.0224, 0.142; Tg0 = 2.725; nu0 = 1420.406
# gas temperature: Compton coupling to the CMB with residual ionisation x_e(z) (simple freeze-out fit), integrated from z = 1000
H0 = 100*0.677*1e5/3.086e24
def Hz(z): return H0*np.sqrt(wm/0.677**2*(1 + z)**3 + 9e-5*(1 + z)**4)
def xe(z): return 2.0e-4 + 0.1*np.exp(-(1000 - min(z, 1000))/60)*(z > 800)
sT, aR, me_c = 6.652e-25, 7.566e-15, 9.109e-28*2.998e10
def dTk(z, T):
    Tg = Tg0*(1 + z); x = xe(z); rate = 8*sT*aR*Tg**4*x/(3*me_c*(1 + 1.08*x))
    return [2*T[0]/(1 + z) - rate*(Tg - T[0])/(Hz(z)*(1 + z))]
_z = np.linspace(1000, 15, 4000); _T = solve_ivp(dTk, [1000, 15], [Tg0*1001], t_eval=_z, rtol=1e-8, max_step=2).y[0]
Tk_of = lambda z: np.interp(z, _z[::-1], _T[::-1])
def Tb(z, d):
    Tg = Tg0*(1 + z); Tk = Tk_of(z)
    nH = 1.9e-7*(1 + z)**3; kap = 3.1e-11*Tk**0.357*np.exp(-32/Tk)
    Ts_ = 0.068*(1 + 2*d); A = 2.85e-15*(1 + 4*d)
    xc = nH*kap*Ts_/(A*Tg); Ts = (Tg + xc*Tk)/(1 + xc)
    return 27*(1 - Tg/Ts)*np.sqrt((1 + z)/10)*(wb/0.023)*np.sqrt(0.15/wm)
print(" z_s   delta   step at (MHz)   T_b there (mK)   jump in T_b (mK)   frequency shift of the feature (kHz)")
for zs in (100, 150, 200):
    for d in (0.004, 0.011):
        nus = nu0/(1 + zs)
        above = Tb(zs, 0.0)                                   # just above nu_s: emitted after the switch
        zb = (1 + zs)*(1 + 2*d) - 1                           # just below nu_s: emitted before, at the shifted line
        below = Tb(zb, d)
        print(f"{zs:5d}  {d:5.3f}   {nus:8.2f}        {above:8.2f}          {below - above:+8.3f}            {2*d*nus*1e3:8.0f}")
print("dark-ages trough for reference: min T_b =", f"{min(Tb(z, 0) for z in np.arange(20, 300)):.1f} mK at z ~ "
      f"{np.arange(20, 300)[np.argmin([Tb(z, 0) for z in np.arange(20, 300)])]}")
