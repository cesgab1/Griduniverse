"""
Bounce instead of inflation (matter/LambdaCDM-bounce family): is the primordial tilt consistent with the measured running?
In a contracting, dust-dominated phase, modes that leave the sound horizon while the effective equation of state is
w_eff = -delta (slightly negative, from dark energy still present) get n_s - 1 = 12 w_eff/(1 + 3 w_eff) ~ -12 delta
(Cai & Wilson-Ewing 2015, eq. 41: n_s = 1 - 12 delta; n_s = 0.965 needs delta ~ 0.003). They note the running is positive but
leave its size open. With a constant dust sound speed, the exit scale is k = |aH|/c_s ~ a^-1/2 in the dust era, and the
dark-energy fraction falls toward the bounce, so delta(k) is fixed once delta(k*) is fixed:
   Lambda:          Omega_DE ~ a^3     -> delta ~ k^-6,    w_DE = -1
   our grid law:    rho_DE ~ |adot|^-1/2 ~ a^(1/4) in the dust era -> Omega_DE ~ a^3.25 -> delta ~ k^-6.5, w_DE = -1 - 1/12
We evaluate n_s(k) and alpha_s = dn_s/dlnk exactly (not only the leading order) on a background with dust + dark energy,
normalised so that n_s(k*) = 0.965, and compare with Planck 2018: alpha_s = -0.0045 +/- 0.0067.
"""
import numpy as np
from scipy.optimize import brentq
def run(law):
    a = np.geomspace(1e-4, 1.0, 200000)                                   # contracting branch, dust-DE equality at a = 1
    if law == "lambda":
        rde = np.ones_like(a); wde = -1.0*np.ones_like(a); H = np.sqrt(a**-3 + rde)
    else:
        D0 = 2**0.25
        a = a[::50]
        H = np.array([brentq(lambda h: h**2 - x**-3 - D0*(x*h)**-0.5, 1e-6, 1e8, xtol=1e-14, rtol=1e-14) for x in a])
        rde = D0*(a*H)**-0.5
        dl = np.gradient(np.log(rde), np.log(a)); wde = -1 - dl/3                       # continuity: d ln rho/d ln a = -3(1+w)
    weff = wde*rde/H**2; ns = 1 + 12*weff/(1 + 3*weff)
    lk = np.log(a*H)                                                      # sound-horizon exit: k ~ |aH| (constant c_s)
    j = np.argmin(abs(ns - 0.965)); al = np.gradient(ns, lk)
    return ns[j], al[j], weff[j], a[j]
for law in ("lambda", "grid"):
    ns, al, w, aj = run(law)
    print(f"{law:7s}: n_s = {ns:.4f} at exit when w_eff = {w:.4f} (DE fraction {abs(w):.4f}); running alpha_s = {al:+.3f} -> "
          f"{(al + 0.0045)/0.0067:.0f} sigma from Planck 2018")
