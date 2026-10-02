"""
Two real orbits, computed with the model's weak-field gravity (which is Newton's law; see paper/figures/CAPTIONS.md, Fig 2-3).
Units: AU, years, GM_sun = 4 pi^2.
(1) nu Octantis: binary A (1.57 Msun) + B (white dwarf, 0.573 Msun), a = 2.609 AU, e = 0.236 (P = 1050 d);
    planet Ab, 2.17 Jupiter masses, a = 1.2445 AU, e = 0.195 around A, inside the binary orbit (arXiv:2609.06885).
    Test: same planet started RETROGRADE vs PROGRADE (coplanar) -> which survives?
(2) Pluto: Sun + Neptune (a 30.07, circular) + Pluto (a 39.48, e 0.249, i 17.1 deg). Pluto crosses Neptune's orbit, yet the two
    never meet: 3:2 resonance. Test: start Pluto at perihelion with Neptune 90 deg away (real configuration) vs with
    Neptune at the same place (no protection) -> closest approach over 50,000 years.
"""
import numpy as np
from scipy.integrate import solve_ivp
G = 4*np.pi**2; MJ = 9.546e-4
def nbody(m, x0, v0, T, n=4000, rtol=1e-10):
    N = len(m); m = np.array(m)
    def f(t, s):
        x = s[:3*N].reshape(N, 3); v = s[3*N:]; a = np.zeros((N, 3))
        for i in range(N):
            d = x - x[i]; r3 = (np.sum(d**2, 1) + 1e-12)**1.5; r3[i] = np.inf
            a[i] = G*np.sum((m/r3)[:, None]*d, 0)
        return np.r_[v, a.ravel()]
    s = solve_ivp(f, [0, T], np.r_[np.ravel(x0), np.ravel(v0)], t_eval=np.linspace(0, T, n), rtol=rtol, atol=1e-12, method="DOP853")
    return s.t, s.y[:3*N].T.reshape(-1, N, 3), s.y[3*N:].T.reshape(-1, N, 3)
def kepler_state(M, a, e, f_peri=True, i_deg=0.0, retro=False):
    r = a*(1 - e); v = np.sqrt(G*M*(1 + e)/r); x = np.array([r, 0, 0]); u = np.array([0, -v if retro else v, 0])
    ci, si = np.cos(np.radians(i_deg)), np.sin(np.radians(i_deg)); R = np.array([[1, 0, 0], [0, ci, -si], [0, si, ci]])
    return R @ x, R @ u

print("(1) nu Octantis: planet inside a tight binary")
mA, mB, mp = 1.57, 0.573, 2.17*MJ
for retro in (True, False):
    # binary: B at apocentre of its orbit relative to A (planet starts far from B)
    ab, eb = 2.609, 0.236; Mt = mA + mB; ra = ab*(1 + eb); va = np.sqrt(G*Mt*(1 - eb)/ra)
    xB_rel, vB_rel = np.array([-ra, 0, 0]), np.array([0, -va, 0])
    xp_rel, vp_rel = kepler_state(mA + mp, 1.2445, 0.195, retro=retro)
    m = [mA, mB, mp]; x = np.array([[0, 0, 0], xB_rel, xp_rel]); v = np.array([[0, 0, 0], vB_rel, vp_rel])
    cm = (np.array(m)[:, None]*x).sum(0)/sum(m); cv = (np.array(m)[:, None]*v).sum(0)/sum(m); x -= cm; v -= cv
    t, X, V = nbody(m, x, v, 300.0, n=30000)
    d = np.linalg.norm(X[:, 2] - X[:, 0], axis=1); bad = np.where((d > 2.0) | (d < 0.05))[0]
    print(f"   {'RETROGRADE' if retro else 'PROGRADE  '}: distance from star A {d.min():.2f}-{d.max():.2f} AU over 300 yr (~104 binary orbits); "
          + ("STAYS BOUND" if len(bad) == 0 else f"LOST after {t[bad[0]]:.1f} yr"))

print("\n(2) Pluto and Neptune")
for lab, phase in (("real (Neptune 90 deg ahead when Pluto is at perihelion)", np.radians(90)), ("unprotected (Neptune at Pluto's perihelion)", 0.0)):
    aN = 30.07; vN = np.sqrt(G/aN); xN = aN*np.array([np.cos(phase), np.sin(phase), 0]); uN = vN*np.array([-np.sin(phase), np.cos(phase), 0])
    xP, uP = kepler_state(1.0, 39.48, 0.249, i_deg=17.1)
    # put Pluto's perihelion near Neptune's orbit radius direction x (perihelion at 29.6 AU, inside Neptune's 30.07 AU)
    m = [1.0, 5.15e-5, 6.6e-9]; x = np.array([[0, 0, 0], xN, xP]); v = np.array([[0, 0, 0], uN, uP])
    cm = (np.array(m)[:, None]*x).sum(0)/sum(m); cv = (np.array(m)[:, None]*v).sum(0)/sum(m); x -= cm; v -= cv
    t, X, V = nbody(m, x, v, 50000.0, n=60000, rtol=1e-9)
    dNP = np.linalg.norm(X[:, 2] - X[:, 1], axis=1)
    rP = np.linalg.norm(X[:, 2] - X[:, 0], axis=1)
    print(f"   {lab}: closest Pluto-Neptune approach over 50,000 yr = {dNP.min():.1f} AU; Pluto's distance from Sun {rP.min():.1f}-{rP.max():.1f} AU")
