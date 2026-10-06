"""Iteration 114: emergent gravity (Hossenfelder static limit; Verlinde local variant) vs the Bullet Cluster (PREREG_114.md)."""
import os, sys, time, numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, "../bullet_tip/iter102_tip.py")).read().split("# ---------- A: Newton")[0]
src = src.replace('here = os.path.dirname(os.path.abspath(__file__))', '')
exec(src)
Gk = 4.30e-6; CONV = 3.086e19/1e6                          # (km/s)^2/kpc per m/s^2
aV = 2.998e8*(70e3/3.086e22)/6*CONV; a0 = 1.2e-10*CONV
print(f"a_V = {aV/CONV:.3e} m/s^2", flush=True)
emain, esub = np.hypot(.06, .06), np.hypot(.05, .06)
dM = pts["main BCG"][6] - pts["main plasma"][6]; dS = pts["sub BCG"][6] - pts["sub plasma"][6]
Mtot = sum(M for _, M, _ in blobs); cen = sum(c_*M for c_, M, _ in blobs)/Mtot
def fields(gx, gy, gz):
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij"); g = [np.zeros_like(X) for _ in range(3)]; rho = np.zeros_like(X)
    for c_, M, a in blobs:
        dx, dy, dz = X - c_[0], Y - c_[1], Z; s2 = dx*dx + dy*dy + dz*dz + a*a; f = -Gk*M/s2**1.5
        g[0] += f*dx; g[1] += f*dy; g[2] += f*dz; rho += 3*M/(4*np.pi*a**3)*(1 + (dx*dx + dy*dy + dz*dz)/a**2)**-2.5
    return X, Y, Z, g, rho
def score(label, gx, gy, gz, h, gD=None, Sig_ph=None):
    if Sig_ph is None:
        div = sum(np.gradient(gD[i], h, axis=i) for i in range(3)); Sig_ph = (-div/(4*np.pi*Gk)).sum(axis=2)*h
    ip = RegularGridInterpolator((gx, gy), Sig_ph, bounds_error=False, fill_value=0)
    kph = lambda p: np.sum(ip(np.c_[p[0] + ax, p[1] + ay])*aw)/(np.pi*100**2)/Scrit
    k = {n: kap_ap(P[n]) + kph(P[n]) for n in names}
    XB, YB = np.meshgrid(gx, gy, indexing="ij")
    tot = sum(sig_plummer(XB - c_[0], YB - c_[1], M, a) for c_, M, a in blobs) + Sig_ph
    i, j = np.unravel_index(np.argmax(tot), tot.shape); tip = np.array([gx[i], gy[j]])
    near = min(names, key=lambda n: np.linalg.norm(P[n] - tip))
    m, s = k["main BCG"] - k["main plasma"], k["sub BCG"] - k["sub plasma"]
    zm, zs = (dM - m)/emain, (dS - s)/esub
    lines = [f"{label}:", "   kappa: " + ", ".join(f"{n} {k[n]:.3f} (meas {pts[n][6]:.2f})" for n in names),
             f"   galaxy-minus-gas: main {m:+.3f} (meas {dM:+.2f}) -> {zm:.1f} sigma; sub {s:+.3f} (meas {dS:+.2f}) -> {zs:.1f} sigma; combined {np.hypot(zm, zs):.1f} sigma",
             f"   TIP at ({tip[0]:.0f},{tip[1]:.0f}) kpc -> nearest {near} ({np.linalg.norm(P[near]-tip):.0f} kpc)"]
    print("\n".join(lines), flush=True); return lines, k
out = ["Iteration 114 -- emergent gravity vs the Bullet Cluster (visible matter of iteration 102)", ""]
# ---------- fine grid (as iteration 102) ----------
h = 20.0; gx = np.arange(-700, 1201, h); gy = np.arange(-600, 601, h); gz = np.arange(-3000, 3001, h)
X, Y, Z, g, rho = fields(gx, gy, gz); gm = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2)
# A1 curl-free Hossenfelder: |gD| gD = aV gN  -> gD = sqrt(aV/|gN|) gN
fac = np.sqrt(aV/gm); out += score("A1 Hossenfelder static limit, curl-free (20 kpc grid)", gx, gy, gz, h, [fac*g[i] for i in range(3)])[0]
# B Verlinde local: |gD|^2 = aV (|gN| + 4 pi G rho r), along gN
rr = np.sqrt((X - cen[0])**2 + (Y - cen[1])**2 + Z**2)
facB = np.sqrt(aV*(gm + 4*np.pi*Gk*rho*rr))/gm
out += score(f"B Verlinde spread term made local, centre ({cen[0]:.0f},{cen[1]:.0f}) (20 kpc grid)", gx, gy, gz, h, [facB*g[i] for i in range(3)])[0]
# reference: galaxy rule (iteration 102) on the same grid
nu1 = 1/(1 - np.exp(-np.sqrt(gm/a0))) - 1
out += score("reference: galaxy rule QUMOND (iteration 102 method)", gx, gy, gz, h, [nu1*g[i] for i in range(3)])[0]
del X, Y, Z, g, rho, gm, rr, fac, facB, nu1
# ---------- A2: full non-linear solution on a coarse grid ----------
hc = 40.0; cx = np.arange(-1200, 1901, hc); cy = np.arange(-1200, 1201, hc); cz = np.arange(-1200, 1201, hc)
X, Y, Z, g, rho = fields(cx, cy, cz); gm = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2)
fac = np.sqrt(aV/gm); lA1, kA1 = score("A1 curl-free on the coarse box (40 kpc, |z| < 1.2 Mpc) -- baseline for the curl check", cx, cy, cz, hc, [fac*g[i] for i in range(3)])
rb = np.sqrt((X - cen[0])**2 + (Y - cen[1])**2 + Z**2) + 1e-9
phi0 = np.sqrt(Gk*Mtot*aV)*np.log(rb)                    # deep-limit point-mass potential: boundary + start
inner = (slice(1, -1),)*3; shape = phi0.shape; src4 = 4*np.pi*Gk*aV*rho; scale = np.sqrt(Gk*Mtot*aV)
def unpack(v): p = phi0.copy(); p[inner] = v.reshape(tuple(s - 2 for s in shape))*scale; return p
def energy(v):
    p = unpack(v); gxx = np.diff(p, axis=0)[:, :-1, :-1]/hc; gyy = np.diff(p, axis=1)[:-1, :, :-1]/hc; gzz = np.diff(p, axis=2)[:-1, :-1, :]/hc
    m = np.sqrt(gxx**2 + gyy**2 + gzz**2 + 1e-12)
    E = np.sum(m**3)/3 + np.sum(src4*p)
    G_ = np.zeros_like(p); wx, wy, wz = m*gxx/hc, m*gyy/hc, m*gzz/hc
    G_[:-1, :-1, :-1] -= wx + wy + wz; G_[1:, :-1, :-1] += wx; G_[:-1, 1:, :-1] += wy; G_[:-1, :-1, 1:] += wz
    G_ += src4
    return E/scale**3*hc**3, G_[inner].ravel()*scale/scale**3*hc**3
t = time.time()
res = minimize(energy, phi0[inner].ravel()/scale, jac=True, method="L-BFGS-B", options=dict(maxiter=6000, maxcor=30, gtol=1e-10, ftol=1e-15))
print(f"A2 solve: {res.nit} iterations, {time.time()-t:.0f} s, {res.message}", flush=True)
p = unpack(res.x)
lap = sum(np.gradient(np.gradient(p, hc, axis=i), hc, axis=i) for i in range(3))
Sig = (lap/(4*np.pi*Gk))[:, :, 1:-1].sum(axis=2)*hc
lA2, kA2 = score("A2 Hossenfelder static limit, FULL non-linear solution (40 kpc)", cx, cy, cz, hc, Sig_ph=Sig)
kv = {n: kap_ap(P[n]) for n in names}
out += [""] + lA1 + lA2 + ["   curl term: extra kappa A2/A1 at markers: " + ", ".join(f"{n} {(kA2[n]-kv[n])/(kA1[n]-kv[n]):.2f}" for n in names),
        f"   (solver: {res.nit} iterations, {res.message})"]
txt = "\n".join(out); open(os.path.join(here, "iter114_bullet_eg.txt"), "w").write(txt + "\n"); print("\nDONE")
