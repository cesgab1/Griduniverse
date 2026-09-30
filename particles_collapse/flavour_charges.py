"""
A whole-number flavour charge for every particle (Froggatt-Nielsen), filtered by the grid's anomaly rule.
Yukawa couplings:  y_ij = c_ij * eps^n_ij,  c_ij = random O(1) complex numbers (|c| in [0.5, 2]),  eps ~ 0.2
  up:   n = |Q_i + u_j - h|      down: n = |Q_i + d_j + h|      leptons: n = |L_i + e_j + h|
Targets (at M_Z): 9 charged-fermion masses + |V_us|, |V_cb|, |V_ub|. A fit is 'good' if every observable is within a factor ~3.
Test A  continuous U(1) flavour symmetry with only Standard-Model fermions, anomalies cancelled exactly:
        [SU3]^2 U1:  sum(2Q+u+d) = 0     [SU2]^2 U1: sum(3Q+L) = 0
        [Y]^2 U1:    sum(Q+8u+2d+3L+6e) = 0      Y [U1]^2: sum(Q^2-2u^2+d^2-L^2+e^2) = 0      (U1^3, gravity: fixed by singlets)
Test B  discrete Z_N flavour symmetry (Z_N link variables, the simplest gauge field a grid can carry):
        exponents n -> min(n mod N, N - n mod N);  robust anomaly rule only: sum(2Q+u+d) = 0 and sum(3Q+L) = 0  (mod N)
"""
import numpy as np, sys, json
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
v = 246.22
m = dict(u=1.29e-3, c=0.619, t=171.7, d=2.9e-3, s=0.055, b=2.89, e=0.4866e-3, mu=0.1027, tau=1.746)   # running masses at M_Z (GeV)
obs = np.log(np.array([m["u"], m["c"], m["t"], m["d"], m["s"], m["b"], m["e"], m["mu"], m["tau"]]) * np.sqrt(2) / v)
obs = np.r_[obs, np.log([0.2250, 0.0418, 0.00369])]
obs_nu = np.log([0.307, 0.52, 0.0220, 0.0296])
ND = 48
import os
NU = int(os.environ.get('NU', 0))
def draws():
    mag = np.exp(rng.uniform(np.log(0.5), np.log(2), (4, ND, 3, 3))); ph = np.exp(2j*np.pi*rng.random((4, ND, 3, 3)))
    return mag * ph
C = draws()

def predict(ch, eps, N=None):
    Q, u, d, L, e, h = ch[0:3], ch[3:6], ch[6:9], ch[9:12], ch[12:15], ch[15]
    nu = np.abs(Q[:, None] + u[None, :] - h); nd = np.abs(Q[:, None] + d[None, :] + h); ne = np.abs(L[:, None] + e[None, :] + h)
    if N:
        f = lambda n: np.minimum(n % N, N - n % N)
        nu, nd, ne = f(Q[:, None] + u[None, :] - h), f(Q[:, None] + d[None, :] + h), f(L[:, None] + e[None, :] + h)
    Yu, Yd, Ye = C[0] * eps**nu, C[1] * eps**nd, C[2] * eps**ne
    Uu, su, _ = np.linalg.svd(Yu); Ud, sd, _ = np.linalg.svd(Yd); se = np.linalg.svd(Ye, compute_uv=False)
    V = np.abs(np.einsum("kji,kjl->kil", Uu.conj(), Ud))          # U_u^dagger U_d ; svd sorts descending -> index 0 = heaviest
    out = np.c_[np.log(su[:, ::-1]), np.log(sd[:, ::-1]), np.log(se[:, ::-1]), np.log(V[:, 2, 1]), np.log(V[:, 1, 0]), np.log(V[:, 2, 0])]
    # V index: heaviest=0 -> (u,c,t)=(2,1,0); V_us = V[c-row? ] use generation order: V_us = |V|[u=2, s=1], V_cb = [c=1, b=0], V_ub = [u=2, b=0]
    if NU:
        f = (lambda n: np.minimum(n % N, N - n % N)) if N else np.abs
        nn = f(L[:, None] + L[None, :] + 2*h)
        Cn = (C[3] + np.transpose(C[3], (0, 2, 1))) / 2
        Un, sn, _ = np.linalg.svd(Cn * eps**nn); Ue, _, _ = np.linalg.svd(Ye)
        Ue = Ue[:, :, ::-1]; Un = Un[:, :, ::-1]; sn = sn[:, ::-1]
        U = np.abs(np.einsum("kji,kjl->kil", Ue.conj(), Un))**2
        s13 = U[:, 0, 2]; s12 = U[:, 0, 1] / (1 - s13 + 1e-12); s23 = U[:, 1, 2] / (1 - s13 + 1e-12)
        r = (sn[:, 1]**2 - sn[:, 0]**2) / (sn[:, 2]**2 - sn[:, 0]**2 + 1e-300)
        mlight = sn[:, 0] / sn[:, 2]
        Up = np.einsum('kji,kjl->kil', Ue.conj(), Un); mbb = np.abs(np.einsum('ki,ki->k', Up[:, 0, :]**2, sn)) / sn[:, 2]
        out = np.c_[out, np.log(s12 + 1e-12), np.log(s23 + 1e-12), np.log(s13 + 1e-12), np.log(np.abs(r) + 1e-12), np.log(mlight + 1e-30), np.log(mbb + 1e-30)]
    return np.median(out, axis=0)

def fixV(p):   # reorder the CKM entries computed above into (V_us, V_cb, V_ub)
    return p
def score(ch, eps, N=None):
    p = predict(ch, eps, N)
    sc = np.sum(((p[:12] - obs) / np.log(3.0))**2)
    if NU: sc += np.sum(((p[12:16] - obs_nu) / np.log(1.6))**2)     # neutrino angles are O(1): demand ~60% accuracy
    return sc, p

def anomalies(ch):
    Q, u, d, L, e = ch[0:3], ch[3:6], ch[6:9], ch[9:12], ch[12:15]
    return np.array([np.sum(2*Q + u + d), np.sum(3*Q + L), np.sum(Q + 8*u + 2*d + 3*L + 6*e), np.sum(Q**2 - 2*u**2 + d**2 - L**2 + e**2)])

EPSMIN = float(__import__("os").environ.get("EPSMIN", 0.1))
def anneal(mode, N=None, steps=int(__import__("os").environ.get("STEPS", 40000)), R=int(__import__("os").environ.get("RANGE", 6))):
    ch = rng.integers(-3, 4, 16); eps = 0.22
    def cost(c, ep):
        s, _ = score(c, ep, N)
        a = anomalies(c)
        if mode == "U1": pen = 25 * np.sum(np.abs(a))
        elif mode == "ZN": pen = 25 * (int(a[0] % N != 0) + int(a[1] % N != 0))
        else: pen = 0
        return s + pen
    cur = cost(ch, eps); best = (cur, ch.copy(), eps); T = 20.0
    for it in range(steps):
        T = 20.0 * (0.02 / 20.0)**(it / steps)
        new = ch.copy(); k = rng.integers(16); new[k] = np.clip(new[k] + rng.choice([-1, 1]), -R, R)
        if rng.random() < 0.3:          # also move a compensating charge to walk along anomaly-free directions
            k2 = rng.integers(16); new[k2] = np.clip(new[k2] + rng.choice([-1, 1]), -R, R)
        ne = float(np.clip(eps + rng.normal(0, 0.01), EPSMIN, 0.35))
        c = cost(new, ne)
        if c < cur or rng.random() < np.exp(-(c - cur) / T):
            ch, eps, cur = new, ne, c
            if c < best[0]: best = (c, ch.copy(), eps)
    return best

names = ["u", "c", "t", "d", "s", "b", "e", "mu", "tau", "V_us", "V_cb", "V_ub"]
def report(tag, best, N=None):
    c, ch, eps = best; s, p = score(ch, eps, N); a = anomalies(ch)
    worst = np.max(np.abs(p[:12] - obs)); ratio = np.exp(p[:12] - obs)
    if NU: print(f"   neutrinos: sin2th12 {np.exp(p[12]):.3f} (0.307) sin2th23 {np.exp(p[13]):.3f} (0.52) sin2th13 {np.exp(p[14]):.4f} (0.022) dm21/dm31 {np.exp(p[15]):.3f} (0.030) | PREDICTED m1/m3 = {np.exp(p[16]):.2e}, m_bb/m3 = {np.exp(p[17]):.3f}")
    print(f"{tag}: fit score {s:.2f} (0 = perfect; every observable within x3 needs roughly < 12), worst factor {np.exp(worst):.1f}, eps = {eps:.3f}")
    print(f"   charges Q={ch[0:3]} u={ch[3:6]} d={ch[6:9]} L={ch[9:12]} e={ch[12:15]} h={ch[15]} | anomalies [SU3^2, SU2^2, Y^2, Y*F^2] = {a}")
    print("   predicted / measured: " + ", ".join(f"{n} {r:.2f}" for n, r in zip(names, ratio)))
    return dict(score=float(s), eps=float(eps), charges=ch.tolist(), anomalies=a.tolist(), ratios=ratio.tolist(), worst=float(np.exp(worst)))

mode = sys.argv[2] if len(sys.argv) > 2 else "free"
N = int(sys.argv[3]) if len(sys.argv) > 3 else None
bests = [anneal(mode, N) for _ in range(4)]
best = min(bests, key=lambda b: b[0])
r = report(f"{mode}{'' if N is None else ' N='+str(N)}", best, N)
json.dump(r, open(f"/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/flav_{mode}_{N}_{__import__('os').environ.get('TAG','')}.json", "w"))
