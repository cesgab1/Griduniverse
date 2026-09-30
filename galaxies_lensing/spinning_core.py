"""
Spinning capped black hole (rotating Hayward core, Newman-Janis construction, Bambi & Modesto 2013) vs Kerr.
Kerr with the mass replaced by m(r) = M r^3 / (r^3 + 2 M l^2): the grid's cap l makes the core finite.  Units G = c = M = 1.
Computed on the equator: horizon, horizon spin rate Omega_H (sets jet power, Blandford-Znajek ~ Omega_H^2),
innermost stable circular orbit (ISCO), energy released by matter spiralling in (efficiency = 1 - E_isco),
prograde photon orbit, and frame-dragging 'swirl' speed of the grid.
Framework value of l for real black holes: l ~ 0.6 Planck lengths x kappa, i.e. l/M ~ 1e-38 (stellar) to 1e-48 (M87*).
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
def m(r, l): return r**3/(r**3 + 2*l**2) if l > 0 else 1.0 + 0*r
def metric(r, a, l):
    mm = m(r, l)
    gtt = -(1 - 2*mm/r); gtp = -2*mm*a/r; gpp = r**2 + a**2 + 2*mm*a**2/r
    return gtt, gtp, gpp
def d(f, r, h=1e-6): return (f(r+h) - f(r-h))/(2*h)
def horizon(a, l):
    D = lambda r: r**2 - 2*m(r, l)*r + a**2
    rs = np.linspace(0.05, 3, 6000); v = D(rs); idx = np.where(np.diff(np.sign(v)))[0]
    return brentq(D, rs[idx[-1]], rs[idx[-1]+1]) if len(idx) else None
def circ(r, a, l):
    gtt, gtp, gpp = metric(r, a, l)
    dtt = d(lambda x: metric(x, a, l)[0], r); dtp = d(lambda x: metric(x, a, l)[1], r); dpp = d(lambda x: metric(x, a, l)[2], r)
    disc = dtp**2 - dtt*dpp
    if disc < 0: return None
    Om = (-dtp + np.sqrt(disc))/dpp                                   # prograde
    den = -(gtt + 2*gtp*Om + gpp*Om**2)
    if den <= 0: return None
    return -(gtt + gtp*Om)/np.sqrt(den), Om, den
def isco(a, l, rh):
    rs = np.linspace(rh*1.001 + 1e-3, 12, 4000); E = []
    for r in rs:
        c = circ(r, a, l); E.append(c[0] if c else np.nan)
    E = np.array(E); ok = np.isfinite(E) & (E > 0)
    i = np.nanargmin(np.where(ok, E, np.nan)); return rs[i], E[i]
def photon(a, l, rh):
    rs = np.linspace(rh*1.0001, 5, 20000)
    for r in rs:
        c = circ(r, a, l)
        if c is not None: return r
    return np.nan
def omega_H(a, l, rh):
    gtt, gtp, gpp = metric(rh, a, l); return -gtp/gpp
print(f"{'spin a':>7} {'core l':>7} | {'horizon':>8} {'Omega_H':>8} {'ISCO':>6} {'efficiency':>10} {'photon orbit':>12} {'swirl at ISCO':>13}")
res = []
for a in (0.0, 0.5, 0.9, 0.99):
    for l in (0.0, 0.1, 0.3):
        rh = horizon(a, l)
        if rh is None: print(f"{a:7.2f} {l:7.2f} | no horizon (spin too high for this core)"); continue
        ri, Ei = isco(a, l, rh); rp = photon(a, l, rh); OH = omega_H(a, l, rh)
        gtt, gtp, gpp = metric(ri, a, l); sw = -gtp/gpp*ri
        res.append((a, l, rh, OH, ri, 1-Ei, rp))
        print(f"{a:7.2f} {l:7.2f} | {rh:8.3f} {OH:8.4f} {ri:6.2f} {100*(1-Ei):9.1f}% {rp:12.3f} {sw:12.3f}c")
# maximum spin allowed with a horizon
print("\nmaximum spin that still has a horizon:")
for l in (0.0, 0.1, 0.3, 0.5):
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo+hi)/2
        (lo, hi) = (mid, hi) if horizon(mid, l) else (lo, mid)
    print(f"   core l = {l:.1f}: a_max = {lo:.4f}")
# scaling to real black holes
a = 0.9; base = [r for r in res if r[0] == a and r[1] == 0][0]; small = [r for r in res if r[0] == a and r[1] == 0.1][0]
frac = abs(small[5]/base[5] - 1)
print(f"\nat spin 0.9 a core of l = 0.1 changes the efficiency by {frac*100:.2f}% (effect scales as l^2)")
for name, lM in (("stellar black hole (10 Msun)", 0.6*1.6e-35/1.48e4), ("Sgr A* (4e6 Msun)", 0.6*1.6e-35/5.9e9), ("M87* (6.5e9 Msun)", 0.6*1.6e-35/9.6e12)):
    print(f"   {name}: l/M = {lM:.0e} -> fractional change ~ {frac*(lM/0.1)**2:.0e}")
