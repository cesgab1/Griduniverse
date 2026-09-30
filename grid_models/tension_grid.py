"""
Gravity as the tension of the grid.
Random 3D grid (Poisson points, density 1, Voronoi-neighbour links), every link a string with tension tau.
A mass pushes on the centre point with force F in a 4th direction; each point is displaced by h (the 'dimple').
Link energy: tau (h_i - h_j)^2 / (2 l_ij). Optional local anchoring k h_i^2 / 2 at every point (a 'loose' grid: held in place locally,
little long-range tension). Equilibrium:   sum_j (tau/l_ij)(h_i - h_j) + k h_i = F delta_i,centre,  h = 0 on the outer sphere.
Continuum prediction: T_eff = tau * sum_links l / (3V)  (effective sheet tension), and
   k = 0 :  h(r) = F/(4 pi T_eff) (1/r - 1/R)        -> Newton's 1/r law, with 'G' = 1/(4 pi T_eff): more tension = weaker gravity
   k > 0 :  h(r) ~ exp(-r/lambda)/r, lambda = sqrt(T_eff/k)  -> the pull only reaches a distance set by tension vs anchoring
Waves: point masses mu on the nodes -> speed v = sqrt(T_eff / (mu n)).
"""
import numpy as np, json
from scipy.spatial import Voronoi
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
rng = np.random.default_rng(3)
Rb, pad = 12.0, 3.0
Lh = Rb + pad; n = rng.poisson((2*Lh)**3); P = rng.uniform(-Lh, Lh, (n, 3))
vor = Voronoi(P)
I, J = vor.ridge_points[:, 0], vor.ridge_points[:, 1]
ell = np.linalg.norm(P[I] - P[J], axis=1)
c = np.argmin(np.linalg.norm(P, axis=1)); X = P - P[c]; r = np.linalg.norm(X, axis=1)
inside = r < Rb; idx = -np.ones(n, int); idx[inside] = np.arange(inside.sum()); N = inside.sum()
both = inside[I] & inside[J]
T_eff_per_tau = ell[both].sum() / (3 * (4/3) * np.pi * Rb**3) + 0.5 * ell[inside[I] ^ inside[J]].sum() / (3 * (4/3) * np.pi * Rb**3)
def solve(tau, k, F=1.0):
    w = tau / ell
    rows, cols, vals = [], [], []; diag = np.full(N, k, float)
    for a, b, ww in zip(I, J, w):
        for s, t in ((a, b), (b, a)):
            if inside[s]:
                diag[idx[s]] += ww
                if inside[t]: rows.append(idx[s]); cols.append(idx[t]); vals.append(-ww)
    A = coo_matrix((np.r_[vals, diag], (np.r_[rows, np.arange(N)], np.r_[cols, np.arange(N)])), (N, N)).tocsr()
    f = np.zeros(N); f[idx[c]] = F
    return spsolve(A, f)
rr = r[inside]; sel = (rr > 2.5) & (rr < Rb - 1.5)
print(f"grid: {N} points inside radius {Rb}, effective sheet tension T_eff = {T_eff_per_tau:.3f} x link tension")
print("\n1. Tension only (no anchoring): is the dimple Newton's 1/r, with strength 1/(4 pi T_eff)?")
res = {}
for tau in (0.5, 1.0, 2.0, 4.0):
    h = solve(tau, 0.0)
    Xf = np.c_[1/rr[sel], np.ones(sel.sum())]; coef, b = np.linalg.lstsq(Xf, h[sel], rcond=None)[0]
    pwr = -np.polyfit(np.log(rr[sel]), np.log(h[sel] - b), 1)[0]
    pred = 1 / (4*np.pi*T_eff_per_tau*tau)
    print(f"   link tension {tau:3.1f}: dimple = {coef:.4f}/r (predicted {pred:.4f}/r, ratio {coef/pred:.3f}); falls off as r^-{pwr:.2f}")
    res[f"tau{tau}"] = dict(coef=coef, pred=pred, power=pwr)
print("\n2. Loose grid: tension weakened relative to local anchoring (k = 0.05): how far does the pull reach?")
for tau in (4.0, 1.0, 0.25, 0.05):
    h = solve(tau, 0.05)
    bins = np.arange(2.5, Rb - 1.5, 1.0); prof = [(lo + .5, np.mean(h[(rr >= lo) & (rr < lo+1)])) for lo in bins]
    x = np.array([p[0] for p in prof]); y = np.array([p[1] for p in prof])
    ok = y > 1e-12
    lam = -1 / np.polyfit(x[ok], np.log(y[ok] * x[ok]), 1)[0] if ok.sum() > 3 else float("nan")
    print(f"   link tension {tau:4.2f}: pull range lambda = {lam:5.2f} grid spacings (predicted sqrt(T_eff/k) = {np.sqrt(T_eff_per_tau*tau/0.05):5.2f});"
          f" dimple at r=3: {np.interp(3, x, y):.2e}, at r=9: {np.interp(9, x, y):.2e}")
print("\n3. Waves: speed of ripples through the stretched grid (node mass mu = 1)")
inb = r < Rb
for tau in (1.0, 4.0):
    vs = []
    for _ in range(8):
        kv = rng.normal(size=3); kv *= 0.15 / np.linalg.norm(kv)
        d = X[J] - X[I]; ok = inb[I] & inb[J]
        num = np.sum((tau/ell[ok]) * 2 * (1 - np.cos(d[ok] @ kv))); den = inb.sum()
        vs.append(np.sqrt(num / den) / 0.15)
    print(f"   link tension {tau}: wave speed {np.mean(vs):.3f} (predicted sqrt(T_eff/rho) = {np.sqrt(T_eff_per_tau*tau/1.0):.3f})")
G, cl, lP = 6.6743e-11, 2.99792458e8, 1.616255e-35
T_needed = cl**4 / (4*np.pi*G); ld = 0.603*lP
tau_link = T_needed / T_eff_per_tau * ld**2 / ld**2      # T_eff (N) = tau * [sum l/(3V)] with lengths in units of l_d -> tau in N
print(f"\n4. Real numbers: Newton's G requires sheet tension T = c^4/(4 pi G) = {T_needed:.2e} N; with the grid's link density that is"
      f" {T_needed/T_eff_per_tau:.2e} N per link = {T_needed/T_eff_per_tau/(cl**4/G):.3f} Planck forces (c^4/G)")
print(f"   Dark energy's tension (negative pressure) today: {5.36e-10:.1e} Pa, i.e. {5.36e-10*ld**2:.1e} N per grid-cell face: {5.36e-10*ld**2/(T_needed/T_eff_per_tau):.0e} of the link tension")
json.dump(res, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/tension_results.json", "w"))
