"""
ITERATION 52 (pre-registered in PREREG_52.md): coefficients a, b of the eps^2 (a R^2 + b Ric^2) corrections of our curvature
reader, from exact optimal transport between finely discretised balls on two homogeneous spaces.
"""
import numpy as np, ot, sys
# REVISION after the first run (recorded in iter52_first_runs.txt): a hard-edged lattice ball has the wrong second moment, which
# rescales the leading term by a resolution-dependent constant (kappa/eps^2 -> 0.198 at m = 6, 0.196 at m = 9) and the eps-fit
# then turned that constant into spurious eps^2/eps^3 terms. Fixes: (1) cells weighted by the fraction of their volume inside the
# ball (sub-sampled), (2) the fit includes a free overall factor f (discretisation of a fixed pattern rescales the leading term
# by an eps-independent amount), (3) results compared across resolutions m.
def ball_pattern(eps, m, sub=6):
    g = np.arange(-m - 1, m + 2)/m*eps; X = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    X = X[np.linalg.norm(X, axis=1) < eps + np.sqrt(3)*eps/m]
    o = (np.arange(sub) + 0.5)/sub - 0.5; O = np.stack(np.meshgrid(o, o, o, indexing="ij"), -1).reshape(-1, 3)*eps/m
    frac = np.array([np.mean(np.linalg.norm(x + O, axis=1) < eps) for x in X]); keep = frac > 0
    return X[keep], frac[keep]
def kappa(space, axis, eps, m):
    t, fr = ball_pattern(eps, m); r = np.linalg.norm(t, axis=1); d = eps/2
    if space == "S3":
        w = np.sinc(r/np.pi)**2
        dirn = np.zeros((len(t), 3)); nz = r > 0; dirn[nz] = t[nz]/r[nz, None]
        P = np.column_stack([np.cos(r), np.sin(r)[:, None]*dirn])
        Rot = np.eye(4); Rot[0, 0] = Rot[1, 1] = np.cos(d); Rot[1, 0] = np.sin(d); Rot[0, 1] = -np.sin(d)
        Q = P @ Rot.T; C = np.arccos(np.clip(P @ Q.T, -1, 1))
    else:   # S2 x line: coordinates (u1, u2) on the sphere, w on the line
        u = t[:, :2]; ru = np.linalg.norm(u, axis=1); w = np.sinc(ru/np.pi)
        dirn = np.zeros((len(t), 2)); nz = ru > 0; dirn[nz] = u[nz]/ru[nz, None]
        S = np.column_stack([np.sin(ru)[:, None]*dirn, np.cos(ru)]); Z = t[:, 2]
        if axis == "sphere":
            Ry = np.array([[np.cos(d), 0, np.sin(d)], [0, 1, 0], [-np.sin(d), 0, np.cos(d)]]); S2, Z2 = S @ Ry.T, Z
        else:
            S2, Z2 = S, Z + d
        C = np.sqrt(np.arccos(np.clip(S @ S2.T, -1, 1))**2 + (Z[:, None] - Z2[None, :])**2)
    w = w*fr; w = w/w.sum()
    return 1 - ot.emd2(w, w, C, numItermax=10_000_000)/d
def flat_kappa(eps, m):
    t, fr = ball_pattern(eps, m); w = fr/fr.sum(); d = eps/2
    C = np.linalg.norm(t[:, None, :] - (t + [d, 0, 0])[None, :, :], axis=2); return 1 - ot.emd2(w, w, C, numItermax=10_000_000)/d
m = int(sys.argv[1]) if len(sys.argv) > 1 else 9
out = [f"ITERATION 52: curvature-squared coefficients (expectations pre-registered); lattice m = {m} per ball radius", ""]
out.append(f"C1 flat check: kappa = {flat_kappa(0.3, m):+.2e}")
epss = np.array([0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6])
rows = []
for e in epss:
    k3 = kappa("S3", None, e, m); ks = kappa("S2R", "sphere", e, m); kl = kappa("S2R", "line", e, m)
    R3 = 30*k3/e**2; R2 = 30*(2*ks + kl)/3/e**2; rows.append((e, k3, ks, kl, R3, R2))
    out.append(f"eps {e:.2f}: S3 kappa/eps^2 = {k3/e**2:.5f} (lead 0.2); S2xR sphere-axis {ks/e**2:.5f} (0.1), line-axis {kl/e**2:+.5f} (0)"
               f"; R_eps: S3 {R3:.5f} (R=6), S2xR {R2:.5f} (R=2)")
rows = np.array(rows)
def fit(Reps, R):   # R_eps = f R + A eps^2 + B eps^4 (+ C eps^3)
    M = np.column_stack([R*np.ones_like(epss), epss**2, epss**4]); c = np.linalg.lstsq(M, Reps, rcond=None)[0]
    M3 = np.column_stack([R*np.ones_like(epss), epss**2, epss**3, epss**4]); c3 = np.linalg.lstsq(M3, Reps, rcond=None)[0]
    return c, c3
(c3, c3x), (c2, c2x) = fit(rows[:, 4], 6.0), fit(rows[:, 5], 2.0)
for lab, A3, A2, f3, f2 in (("fit f + eps^2 + eps^4", c3[1], c2[1], c3[0], c2[0]), ("with eps^3 too", c3x[1], c2x[1], c3x[0], c2x[0])):
    a, b = np.linalg.solve([[36, 12], [4, 2]], [A3, A2])
    out.append(f"{lab}: f(S3) = {f3:.5f}, f(S2xR) = {f2:.5f}; A(S3) = {A3:+.4f}, A(S2xR) = {A2:+.4f}  ->  a = {a:+.4f}, b = {b:+.4f}")
out.append(f"eps^3 coefficients when allowed: S3 {c3x[2]:+.4f}, S2xR {c2x[2]:+.4f}")
txt = "\n".join(out); print(txt); open(f"iter52_curvature_squared_m{m}.txt", "w").write(txt + "\n")
