"""
3-D particle-mesh simulation of Khronon's dark fluid around a galaxy (Blanchet & Skordis 2024, weak field, one-fluid bookkeeping).
No hand-picked crossover schedule: the fluid law is the DBI form K(Q) with the parameters that passed our CLASS cosmology fits
(1/mu = 300 Mpc, lambda_D = 9.2e6; sigma8 0.790, Lyman-alpha 0.979).

Physical equations (x = Q - 1, Xi = -c^2 x up to a constant, y = x - x_bg):
  Delta phi = 4 pi G (rho_b + rho_tau - mean)                                         (3.11)
  -div((1 - f) grad y) + (1/2)[K_Q(x_bg + y) - K_Q(x_bg)] = (4 pi G / c^2) drho_tau   (3.13, general K)
  fluid acceleration = -grad phi + grad Xi = -grad phi - c^2 grad y                    (Dv/Dt = -grad phi + grad Xi)
  (1/2) K_Q(x) = mu^2 x / sqrt(1 - lambda x^2)  (DBI), background (1/2) K_Q(x_bg) = (4 pi G/c^2) rho_bar
  f = MOND function of |grad Xi|/a0 = c^2 |grad y| / a0  (inverse of nu(y) = 1/(1 - exp(-sqrt y)))
Comoving coordinates, periodic box, CIC mass assignment, FFT Poisson, Newton-CG for y (FFT-preconditioned), KDK leapfrog in ln a.
The fluid is followed with particles (single-stream fluid -> particles, as in every CDM N-body code). Baryons: one fixed Plummer
mass at the box centre (Newtonian gravity only; in Khronon MOND arises from the fluid, so baryons do not source Xi).
Initial conditions (z = 30): Zel'dovich displacements from the conditional-mean profile of a 1.5e12 Msun region collapsing at
z = 1 (peak_env.py) plus a Gaussian random field (BBKS, sigma8 = 0.81), so orbits are not purely radial.
Mode "cdm": the Xi force is switched off (plain dark matter control). Units: kpc (comoving unless noted), Gyr, Msun.
"""
import numpy as np, sys, os, json, time
from numpy.fft import rfftn, irfftn, fftfreq, rfftfreq
kv = 1.0227; G = 4.30e-6*kv**2; cl = 2.998e5*kv; a0 = 1.2e-10*3.086e19/1e6*kv**2
h = 0.68; H0 = 0.1*h*kv; Om, Ob, OL = 0.31, 0.049, 0.69; Oc = Om - Ob
rhoc = 3*H0**2/(8*np.pi*G); rho_f = Oc*rhoc                                   # comoving mean fluid density
MU = 1/3.0e5; LAM = 9.2e6
TAU = float(os.environ.get('TAU', '0'))                                      # optional dissipation time (Gyr), 0 = off
GEXT = float(os.environ.get('GEXT', '0'))                                   # external field in units of a0 (MOND function only)                                                     # DBI: 1/mu = 300 Mpc, lambda_D
nu = lambda y: 1/(1 - np.exp(-np.sqrt(np.maximum(y, 1e-14))))
_ys = np.logspace(-12, 8, 6000); _xs = _ys*nu(_ys); _fs = 1/nu(_ys)
def fmond(x): return np.interp(x, _xs, _fs, left=_fs[0], right=1.0)
Hof = lambda a: H0*np.sqrt(Om/a**3 + OL)

class Box:
    def __init__(s, L, Ng):
        s.L, s.N = L, Ng; s.dx = L/Ng
        k1 = 2*np.pi*fftfreq(Ng, d=s.dx); kz = 2*np.pi*rfftfreq(Ng, d=s.dx)
        s.KX, s.KY, s.KZ = np.meshgrid(k1, k1, kz, indexing="ij")
        # finite-difference Laplacian eigenvalues (consistent with the stencil used for gradients/divergence)
        s.K2 = (2/s.dx*np.sin(s.KX*s.dx/2))**2 + (2/s.dx*np.sin(s.KY*s.dx/2))**2 + (2/s.dx*np.sin(s.KZ*s.dx/2))**2
        s.K2[0, 0, 0] = 1.0
    def cic(s, pos, m):
        N = s.N; g = pos/s.dx - 0.5; i = np.floor(g).astype(int); d = g - i
        rho = np.zeros(N**3)
        for ox in (0, 1):
            wx = d[:, 0] if ox else 1 - d[:, 0]
            for oy in (0, 1):
                wy = d[:, 1] if oy else 1 - d[:, 1]
                for oz in (0, 1):
                    wz = d[:, 2] if oz else 1 - d[:, 2]
                    idx = (((i[:, 0] + ox) % N)*N + (i[:, 1] + oy) % N)*N + (i[:, 2] + oz) % N
                    rho += np.bincount(idx, weights=m*wx*wy*wz, minlength=N**3)
        return rho.reshape(N, N, N)/s.dx**3
    def interp(s, fld, pos):
        N = s.N; g = pos/s.dx - 0.5; i = np.floor(g).astype(int); d = g - i; out = np.zeros(len(pos))
        for ox in (0, 1):
            wx = d[:, 0] if ox else 1 - d[:, 0]
            for oy in (0, 1):
                wy = d[:, 1] if oy else 1 - d[:, 1]
                for oz in (0, 1):
                    wz = d[:, 2] if oz else 1 - d[:, 2]
                    out += fld[(i[:, 0] + ox) % N, (i[:, 1] + oy) % N, (i[:, 2] + oz) % N]*wx*wy*wz
        return out
    def grad(s, f):                                                           # central differences
        return [(np.roll(f, -1, ax) - np.roll(f, 1, ax))/(2*s.dx) for ax in range(3)]
    def poisson(s, src):                                                      # solve Lap phi = src (periodic, zero mean)
        fk = rfftn(src); fk /= -s.K2; fk[0, 0, 0] = 0; return irfftn(fk, s=src.shape)
    def divwgrad(s, w, y):                                                    # div(w grad y), face-centred, conservative
        out = np.zeros_like(y)
        for ax in range(3):
            wf = 0.5*(w + np.roll(w, -1, ax)); fl = wf*(np.roll(y, -1, ax) - y)/s.dx
            out += (fl - np.roll(fl, 1, ax))/s.dx
        return out

def bg_state(a):
    """Background: s = (4 pi G/c^2) rho_bar / mu^2 (so (1/2)K_Q(x_bg) = mu^2 s), x_bg = s/sqrt(1 + lam s^2)."""
    return 4*np.pi*G*rho_f/(cl**2*MU**2*a**3)
def y_of_ds(s1, ds):
    """y = x(s1 + ds) - x(s1), x(s) = s/sqrt(1 + lam s^2), without cancellation. Also returns dy/ds at s1 + ds."""
    s2 = s1 + ds; A1 = 1 + LAM*s1**2; A2 = 1 + LAM*s2**2
    y = ds*(s1 + s2)/(np.sqrt(A1)*np.sqrt(A2)*(s2*np.sqrt(A1) + s1*np.sqrt(A2)))
    return y, A2**-1.5
def smooth(B, f, sig_cells=1.0):
    fk = rfftn(f); k2 = B.KX**2 + B.KY**2 + B.KZ**2; return irfftn(fk*np.exp(-0.5*k2*(sig_cells*B.dx)**2), s=f.shape)
def solve_y(B, a, drho_com, y0, tol=1e-3, newton=12, cgmax=200):
    """Solve  -(1/a^2) div(w grad y) + mu^2 ds(y) = (4 pi G/c^2) drho_com / a^3   for y = x - x_bg (DBI: x < x_lim).
    Written as the minimum of a convex energy (w frozen per Picard step):
        E[y] = sum_faces w (dy)^2 / (2 a^2 dx^2) + sum_cells [mu^2 Phi(y) - rhs y],   Phi' = ds(y), Phi'' = 1/y'(y) > 0,
    with a barrier at the DBI limit (B2 = 1 - lam x^2 > 0). Newton steps (SPD Hessian, preconditioned CG) with a line search
    on E: guaranteed descent and no crossing of the limit. All DBI algebra is cancellation-free (B2 = Bb - lam y (2 x_bg + y)).
    Returns y, ds(y) and info. (Replaces the earlier ds-variable Newton, which stalled at late times; y0 is the warm start.)"""
    s1 = bg_state(a); _, yp0 = y_of_ds(s1, 0.0)
    if MU**2/yp0*a**2 > 10*(np.pi/B.dx)**2 or (1/a - 1) > float(os.environ.get("XI_ON_Z", "1e9")):
        return np.zeros_like(y0), np.zeros_like(y0), {"res": 0.0, "it": -1, "skipped": True}
    dr = np.maximum(drho_com, -0.95*rho_f)
    rhs = 4*np.pi*G/cl**2*dr/a**3
    xb = s1/np.sqrt(1 + LAM*s1**2); Bb = 1/(1 + LAM*s1**2); sBb = np.sqrt(Bb)
    def B2(y): return Bb - LAM*y*(2*xb + y)
    def ds_of(y):
        b2 = B2(y); S = np.sqrt(b2) + sBb
        # (xb+y)/sqrt(B2) - xb/sqrt(Bb), cancellation-free
        return (y*sBb + xb*(sBb - np.sqrt(b2)))/(np.sqrt(b2)*sBb)
    def Phi(y):
        b2 = B2(y); S = np.sqrt(b2) + sBb
        return y**2*(xb*LAM*(2*xb + y)/S + sBb)/(S*sBb)
    def wfun(y):
        gy = B.grad(y); gx0 = cl**2*gy[0]/a + GEXT*a0
        gm = np.sqrt(gx0**2 + (cl**2*gy[1]/a)**2 + (cl**2*gy[2]/a)**2); return np.maximum(1 - fmond(gm/a0), 1e-6)
    def energy(y, w):
        if np.any(B2(y) <= 0): return np.inf
        e = 0.0
        for ax in range(3):
            wf = 0.5*(w + np.roll(w, -1, ax)); e += np.sum(wf*(np.roll(y, -1, ax) - y)**2)/(2*a**2*B.dx**2)
        return e + np.sum(MU**2*Phi(y) - rhs*y)
    y = y0.copy()
    if np.any(B2(y) <= 0): y = np.zeros_like(y0)
    if not np.any(y) and MU**2/yp0*a**2 > 0.3*(np.pi/B.dx)**2:                  # cold start in the stiff era: local solution
        loc = np.clip(rhs/MU**2, -s1*(1 - 1e-12), 1e6*s1); y, _ = y_of_ds(s1, loc); y = np.where(B2(y) > 0, y, 0.0)
    scale = np.sqrt(np.mean(rhs**2)) + 1e-300; info = {}
    for it in range(newton):
        w = wfun(y)
        R = -B.divwgrad(w, y)/a**2 + MU**2*ds_of(y) - rhs
        rn = np.sqrt(np.mean(R**2))/scale; info["res"] = float(rn); info["it"] = it
        if rn < tol: break
        b2 = B2(y); D = MU**2*b2**-1.5
        Aop = lambda v: -B.divwgrad(w, v)/a**2 + D*v
        diag = D + 6*w/(a**2*B.dx**2); Dm = np.exp(np.mean(np.log(D)))
        Pk = 1/(B.K2/a**2 + Dm); Sc = np.sqrt((6/(a**2*B.dx**2) + Dm)/diag)
        Pre = lambda r: Sc*irfftn(rfftn(Sc*r)*Pk, s=r.shape)
        v = np.zeros_like(y); r = -R.copy(); z = Pre(r); p = z.copy(); rz = np.vdot(r, z); r0 = np.sqrt(np.vdot(r, r))
        for k in range(cgmax):
            Ap = Aop(p); al = rz/np.vdot(p, Ap); v += al*p; r -= al*Ap
            if np.sqrt(np.vdot(r, r)) < 1e-4*r0: break
            z = Pre(r); rzn = np.vdot(r, z); p = z + (rzn/rz)*p; rz = rzn
        # fraction-to-boundary: no cell may move more than 90% of the way to the DBI limit in one step
        ylim = (Bb/LAM)/(xb + np.sqrt(xb**2 + Bb/LAM))
        v = np.where(v > 0, np.minimum(v, 0.9*(ylim - y)), v)
        E0 = energy(y, w); alpha = 1.0
        for _ls in range(30):
            Et = energy(y + alpha*v, w)
            if Et < E0: break
            alpha *= 0.5
        y = y + alpha*v; info["alpha"] = alpha; info["cg"] = k
    return y, ds_of(y), info

def initial_conditions(L, Np, zi, seed=1, Mhalo=1.5e12, zc=1.0, rand=True):
    import peak_env as PE
    a = 1/(1 + zi); q = (np.arange(Np) + 0.5)*L/Np
    QX, QY, QZ = np.meshgrid(q, q, q, indexing="ij"); Q = np.stack([QX.ravel(), QY.ravel(), QZ.ravel()], 1)
    Bp = Box(L, Np); c0 = L/2
    # linear density at a: conditional-mean profile + random field (both already scaled to z = zi by peak_env's D30 for zi=30)
    rq = np.sqrt((QX - c0)**2 + (QY - c0)**2 + (QZ - c0)**2)
    Rg = np.geomspace(5, L, 400); Dg = PE.Delta_peak(Rg)*(PE.Dgr(zi)/PE.D30)
    d_loc = np.gradient(Rg**3*Dg, Rg)/(3*Rg**2)                              # local density from enclosed mean
    delta = np.interp(rq, Rg, d_loc)
    if rand:
        rng = np.random.default_rng(seed); kk = np.sqrt(Bp.KX**2 + Bp.KY**2 + Bp.KZ**2)/1e-3   # 1/Mpc
        kk[0, 0, 0] = 1.0
        Pk = np.interp(kk, PE.ks, PE.Pk)*(PE.Dgr(zi))**2*1e9                  # (kpc)^3
        wn = rfftn(rng.normal(size=(Np, Np, Np)))
        dk = wn*np.sqrt(Pk/L**3)*Np**1.5; dk[0, 0, 0] = 0
        delta = delta + irfftn(dk, s=(Np, Np, Np))
    delta -= delta.mean()
    # Zel'dovich: psi = -grad(Lap^-1 delta)
    phi = Bp.poisson(delta); gx, gy, gz = Bp.grad(phi)
    psi = -np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)
    fgr = (-1 + np.sqrt(1 + 24*Oc/Om))/4                                     # growth rate with only the fluid clustering
    pos = (Q + psi) % L; u = a**2*Hof(a)*fgr*psi
    m = rho_f*L**3/Np**3
    return pos, u, m, a

def run(mode, L=6000.0, Ng=128, Np=64, zi=30.0, Mb=1e11, eps_b=10.0, nsteps=600, out="pm3d", rand=True, seed=1):
    pos, u, m, a = initial_conditions(L, Np, zi, seed=seed, rand=rand)
    B = Box(L, Ng); c0 = np.array([L/2]*3); y = np.zeros((Ng,)*3); st = {'ds': np.zeros((Ng,)*3)}
    lna = np.linspace(np.log(a), 0.0, nsteps + 1); dl = lna[1] - lna[0]
    marks = {2.0: None, 0.9: None, 0.3: None, 0.0: None}; log = []
    def accel(pos, a, y):
        rho = B.cic(pos, m); drho = smooth(B, rho - rho_f)                    # one-cell smoothing for both phi and Xi
        phi = B.poisson(4*np.pi*G*drho/a)                                     # comoving Poisson
        gphi = B.grad(phi); F = [-g for g in gphi]                           # du/dt = -grad_c phi
        info = {}
        if mode == "khronon":
            y, _dsv, info = solve_y(B, a, drho, st['ds']); st['ds'] = y          # st['ds'] holds the warm-start y
            gy = B.grad(y); F = [F[i] - cl**2*gy[i] for i in range(3)]        # du/dt += -c^2 grad_c y
        acc = np.stack([B.interp(F[i], pos) for i in range(3)], 1)
        d = pos - c0; d -= L*np.round(d/L); rp = a*d                          # baryons: physical Plummer force, times a
        r2 = (rp**2).sum(1) + eps_b**2; acc += a*(-G*Mb*rp/r2[:, None]**1.5)
        return acc, y, info
    t0 = time.time(); k0 = 0
    ck = f"{out}_{mode}_ckpt.npz"
    if os.path.exists(ck):
        d = np.load(ck, allow_pickle=True); pos, u, k0 = d["pos"], d["u"], int(d["k"]); st['ds'] = d["ds"]
        marks.update(d["marks"].item()); log = list(d["log"]); print(f"resumed at step {k0}", flush=True)
    acc, y, info = accel(pos, np.exp(lna[k0]), y)
    for k in range(k0, nsteps):
        a1 = np.exp(lna[k]); am = np.exp(lna[k] + dl/2); a2 = np.exp(lna[k + 1])
        dt1 = (np.log(am) - np.log(a1))/Hof(np.sqrt(a1*am)); dt2 = (np.log(a2) - np.log(am))/Hof(np.sqrt(am*a2))
        u += acc*dt1
        pos = (pos + u/am**2*(dt1 + dt2)) % L
        acc, y, info = accel(pos, a2, y)
        u += acc*dt2
        if TAU > 0 and mode == "khronon" and info and not info.get("skipped", False):
            # optional dissipation test: random (multi-stream) motions relax toward the local mean flow on a time TAU (Gyr),
            # as if the fluid could radiate the energy away ('cooling'); coherent flow is untouched
            mom = [B.cic(pos, m*u[:, i]) for i in range(3)]; rho_ = B.cic(pos, m) + 1e-30
            um = np.stack([B.interp(smooth(B, mom[i])/smooth(B, rho_), pos) for i in range(3)], 1)
            u = um + (u - um)*np.exp(-(dt1 + dt2)/TAU)
        z = 1/a2 - 1
        for zm in marks:
            if marks[zm] is None and z <= zm + 1e-9:
                d = pos - c0; d -= L*np.round(d/L); rph = a2*np.sqrt((d**2).sum(1))
                exc = {R: float(m*np.sum(rph < R) - rho_f/a2**3*4/3*np.pi*R**3) for R in (30, 50, 100, 200, 300, 500, 1000)}
                marks[zm] = exc
                log.append(dict(z=zm, excess=exc, res=info.get("res"), wall=time.time() - t0))
                print(f"[{mode}] z={zm}: " + " ".join(f"{R}:{v:.2e}" for R, v in exc.items()) + f" | Xi res {info.get('res')} | {time.time()-t0:.0f}s", flush=True)
        if os.environ.get("STOPZ") and z <= float(os.environ["STOPZ"]):
            np.savez(f"{out}_{mode}_state_z{os.environ['STOPZ']}.npz", pos=pos, u=u, k=k + 1, ds=st['ds'],
                     marks=np.array(marks, dtype=object), log=np.array(log, dtype=object))
            print("stopped", flush=True); return log
        if k % 25 == 0 and k > k0:
            np.savez(ck, pos=pos, u=u, k=k + 1, ds=st['ds'], marks=np.array(marks, dtype=object), log=np.array(log, dtype=object))
        if k % 10 == 0:
            print(f"[{mode}] step {k}/{nsteps} z={z:.2f} res={info.get('res')} {time.time()-t0:.0f}s", flush=True)
    json.dump(log, open(f"{out}_{mode}.json", "w"), indent=1)
    np.save(f"{out}_{mode}_pos.npy", pos)
    return log

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "cdm"
    kw = dict(nsteps=int(os.environ.get("NSTEPS", 600)), Ng=int(os.environ.get("NG", 128)), Np=int(os.environ.get("NP", 64)))
    run(mode, out=os.environ.get("OUT", "pm3d"), **kw)
