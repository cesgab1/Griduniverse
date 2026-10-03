"""Validation of the line-of-sight machinery: ordinary LCDM late ISW (matter potential decaying under Lambda).
Theta_l(k) = -2 ∫ dt [d j_l(k chi)/dt] Phi(k,t),  Phi = (3/2) Om a^-1 delta_m / k^2,  delta_m = D+(t)/D+(t0) delta_0(k)
Literature: the late ISW adds roughly 10-30% of the large-angle power (D_2 ~ 100-300 muK^2 of ~1000)."""
import numpy as np, camb
from scipy.integrate import quad, cumulative_trapezoid
from scipy.special import spherical_jn
Om, OD, h = 0.315, 0.685, 0.674
pars = camb.set_params(H0=67.4, ombh2=0.0224, omch2=0.120, mnu=0.0, tau=0.054, As=2.1e-9, ns=0.965, lmax=40)
pars.WantTransfer = True; pars.set_matter_power(redshifts=[0.0], kmax=2.0)
res = camb.get_results(pars); khc, _, pk = res.get_matter_power_spectrum(minkh=1e-5, maxkh=1.0, npoints=600)
base = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)["total"][:, 0]
a = np.geomspace(0.05, 1.0, 400); H = np.sqrt(Om*a**-3 + OD); t = cumulative_trapezoid(1/(a*H), a, initial=0); dt = np.gradient(t)
chi = np.array([quad(lambda x: 1/(x*x*np.sqrt(Om*x**-3 + OD)), ai, 1)[0] for ai in a])
Ip = cumulative_trapezoid(1/(a*H)**3, a, initial=0) + quad(lambda x: 1/(x*np.sqrt(Om*x**-3 + OD))**3, 1e-6, a[0])[0]
D = H*Ip; D /= D[-1]
for ell in (2, 3, 5, 10):
    ks = np.geomspace(0.05, 200, 400); th = []
    for k in ks:
        Phi = 1.5*Om*D/a/k**2; th.append(2*np.sum(spherical_jn(ell, k*chi)*np.gradient(Phi, t)*dt))   # direct: 2 ∫ j_l dPhi/dt dt (matter Phi is NOT zero early)
    th = np.array(th); P0 = np.interp(ks/2997.92458/h*h, khc, pk[0])  # k in H0/c -> h/Mpc: k[h/Mpc] = k*H0/c /h ... 
    kh = ks/(2997.92458/h)                                              # 1 (H0/c) = h/2997.9 h/Mpc
    P0 = np.interp(kh, khc, pk[0])*(h/2997.92458)**3                    # (Mpc/h)^3 -> (c/H0)^3
    Cl = 2/np.pi*np.trapezoid(ks**2*th**2*P0, ks); Dl = ell*(ell + 1)*Cl/2/np.pi*(2.7255e6)**2
    print(f"l = {ell:2d}: late-ISW D_l = {Dl:7.1f} muK^2   (total LCDM D_l = {base[ell]:7.1f}; ISW fraction {Dl/base[ell]:.2f})")
