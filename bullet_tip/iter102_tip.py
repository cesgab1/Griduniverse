"""Iteration 102: tip of the Bullet Cluster's curvature from visible matter only (PREREG_102.md)."""
import numpy as np, os
from scipy.optimize import least_squares
from scipy.integrate import quad
here = os.path.dirname(os.path.abspath(__file__))
KPS = 4.413                                     # kpc per arcsec
def radec(ra, dec):
    h, m, s = map(float, ra.split(":")); d, dm, ds = map(float, dec.replace("-", "").split(":"))
    return (h + m/60 + s/3600)*15*3600, -(d + dm/60 + ds/3600)*3600
pts = {"main BCG": ("06:58:35.3", "-55:56:56.3", 5.5, .6, .54, .08, .36, .06),
       "main plasma": ("06:58:30.2", "-55:56:35.9", 6.6, .7, .23, .02, .05, .06),
       "sub BCG": ("06:58:16.0", "-55:56:35.1", 2.7, .3, .58, .09, .20, .05),
       "sub plasma": ("06:58:21.2", "-55:56:30.0", 5.8, .6, .12, .01, .02, .06)}
ra0, de0 = radec(*pts["main BCG"][:2]); cosd = np.cos(np.radians(de0/3600))
P = {k: np.array([-(radec(v[0], v[1])[0] - ra0)*cosd*KPS, (radec(v[0], v[1])[1] - de0)*KPS]) for k, v in pts.items()}
names = list(pts)
# aperture integration points (100 kpc)
rr, th = np.meshgrid(np.linspace(0, 100, 60)[1:] - 100/118, np.linspace(0, 2*np.pi, 90, endpoint=False))
ax, ay = (rr*np.cos(th)).ravel(), (rr*np.sin(th)).ravel(); aw = (rr*(100/59)*(2*np.pi/90)).ravel()
def sig_plummer(dx, dy, M, a): return M*a*a/(np.pi*(a*a + dx*dx + dy*dy)**2)
def ap_mass(blobs, p):
    return sum(np.sum(sig_plummer(p[0] + ax - c[0], p[1] + ay - c[1], M, a)*aw) for c, M, a in blobs)
def fit(centres, col, ecol, abounds):
    def res(q):
        bl = [(P[centres[0]], 10**q[0], q[2]), (P[centres[1]], 10**q[1], q[3])]
        return [(ap_mass(bl, P[n])/1e12 - pts[n][col])/pts[n][ecol] for n in names]
    r = least_squares(res, [12.8, 12.8, sum(abounds)/2, sum(abounds)/2], bounds=([11, 11, abounds[0], abounds[0]], [14.5, 14.5, abounds[1], abounds[1]]))
    return [(P[centres[0]], 10**r.x[0], r.x[2]), (P[centres[1]], 10**r.x[1], r.x[3])], r.fun
gas, rg = fit(["main plasma", "sub plasma"], 2, 3, (30, 800))
st, rs = fit(["main BCG", "sub BCG"], 4, 5, (5, 300))
out = ["Iteration 102 -- Bullet Cluster, visible matter only",
       "positions (kpc from main BCG): " + ", ".join(f"{n} ({P[n][0]:.0f},{P[n][1]:.0f})" for n in names),
       f"BCG-BCG {np.linalg.norm(P['sub BCG']-P['main BCG']):.0f} kpc; main BCG-plasma {np.linalg.norm(P['main plasma']-P['main BCG'])/KPS:.1f} arcsec",
       "gas fit: " + ", ".join(f"M {M:.2e} a {a:.0f}" for _, M, a in gas) + f"  (pulls {np.round(rg, 2)})",
       "star fit: " + ", ".join(f"M {M:.2e} a {a:.0f}" for _, M, a in st) + f"  (pulls {np.round(rs, 2)})"]
# Sigma_crit (flat LCDM Om 0.3, h 0.7; sources z = 1)
c, G, Mpc, Ms = 2.998e8, 6.674e-11, 3.086e22, 1.989e30
def Dc(z): return quad(lambda x: 1/np.sqrt(0.3*(1 + x)**3 + 0.7), 0, z)[0]*c/(70e3/Mpc)
zl, zs = 0.296, 1.0
Dl, Ds, Dls = Dc(zl)/(1 + zl), Dc(zs)/(1 + zs), (Dc(zs) - Dc(zl))/(1 + zs)
Scrit = c**2/(4*np.pi*G)*Ds/(Dl*Dls)/Ms*(Mpc/1e3)**2          # Msun/kpc^2
out.append(f"Sigma_crit = {Scrit:.2e} Msun/kpc^2 (z_s = 1)")
blobs = gas + st
def kap_ap(p, extra=None):
    k = ap_mass(blobs, p)/(np.pi*100**2)/Scrit
    if extra is not None: k += extra(p)
    return k
# ---------- A: Newton ----------
xs = np.arange(-500, 1001, 10.0); ys = np.arange(-500, 501, 10.0); X, Y = np.meshgrid(xs, ys)
SigB = sum(sig_plummer(X - c_[0], Y - c_[1], M, a) for c_, M, a in blobs)
iy, ix = np.unravel_index(np.argmax(SigB), SigB.shape)
tipA = np.array([xs[ix], ys[iy]])
out += ["", "A. Newton/Einstein, visible matter:"]
kA = {n: kap_ap(P[n]) for n in names}
for n in names: out.append(f"   {n:12s} kappa_vis {kA[n]:.3f}   measured {pts[n][6]:.2f} +/- {pts[n][7]:.2f}")
near = min(names, key=lambda n: np.linalg.norm(P[n] - tipA))
out.append(f"   TIP at ({tipA[0]:.0f},{tipA[1]:.0f}) kpc: nearest marker {near} ({np.linalg.norm(P[near]-tipA):.0f} kpc away)")
dmeas = pts["main BCG"][6] - pts["main plasma"][6]; emeas = np.hypot(.06, .06)
out.append(f"   kappa(main BCG) - kappa(main plasma): measured {dmeas:+.2f} +/- {emeas:.3f}; visible {kA['main BCG']-kA['main plasma']:+.3f}"
           f" -> {(dmeas-(kA['main BCG']-kA['main plasma']))/emeas:.1f} sigma")
# ---------- B: galaxy rule (QUMOND), 3-D grid ----------
h = 20.0
gx = np.arange(-700, 1201, h); gy = np.arange(-600, 601, h); gz = np.arange(-3000, 3001, h)
a0 = 1.2e-10*3.086e19/1e6          # (km/s)^2 / kpc
Gk = 4.30e-6
Xg, Yg, Zg = np.meshgrid(gx, gy, gz, indexing="ij")
gN = [np.zeros_like(Xg) for _ in range(3)]
for c_, M, a in blobs:
    dx, dy, dz = Xg - c_[0], Yg - c_[1], Zg
    f = -Gk*M/(dx*dx + dy*dy + dz*dz + a*a)**1.5
    gN[0] += f*dx; gN[1] += f*dy; gN[2] += f*dz
gmag = np.sqrt(gN[0]**2 + gN[1]**2 + gN[2]**2)
nu1 = 1/(1 - np.exp(-np.sqrt(gmag/a0))) - 1
div = sum(np.gradient(nu1*gN[i], h, axis=i) for i in range(3))
rho_ph = -div/(4*np.pi*Gk)
Sig_ph = rho_ph.sum(axis=2)*h
from scipy.interpolate import RegularGridInterpolator
ip = RegularGridInterpolator((gx, gy), Sig_ph, bounds_error=False, fill_value=0)
def kph(p):
    pts_ = np.c_[p[0] + ax, p[1] + ay]
    return np.sum(ip(pts_)*aw)/(np.pi*100**2)/Scrit
SigT = np.array([[0.0]])
# total map on the 2-D grid of B
XB, YB = np.meshgrid(gx, gy, indexing="ij")
SigBB = sum(sig_plummer(XB - c_[0], YB - c_[1], M, a) for c_, M, a in blobs)
tot = SigBB + Sig_ph
i, j = np.unravel_index(np.argmax(tot), tot.shape); tipB = np.array([gx[i], gy[j]])
out += ["", "B. Galaxy rule (QUMOND, a0 = 1.2e-10), visible matter + its phantom curvature:"]
kB = {n: kA[n] + kph(P[n]) for n in names}
for n in names: out.append(f"   {n:12s} kappa {kB[n]:.3f} (phantom part {kph(P[n]):.3f})   measured {pts[n][6]:.2f} +/- {pts[n][7]:.2f}")
near = min(names, key=lambda n: np.linalg.norm(P[n] - tipB))
out.append(f"   TIP at ({tipB[0]:.0f},{tipB[1]:.0f}) kpc: nearest marker {near} ({np.linalg.norm(P[near]-tipB):.0f} kpc away)")
dB = kB["main BCG"] - kB["main plasma"]
out.append(f"   kappa(main BCG) - kappa(main plasma): measured {dmeas:+.2f}; galaxy rule {dB:+.3f} -> {(dmeas-dB)/emeas:.1f} sigma")
dS = (pts["sub BCG"][6] - pts["sub plasma"][6]); eS = np.hypot(.05, .06)
out.append(f"   kappa(sub BCG) - kappa(sub plasma): measured {dS:+.2f}; Newton {kA['sub BCG']-kA['sub plasma']:+.3f}; galaxy rule {kB['sub BCG']-kB['sub plasma']:+.3f}")
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter102_tip.txt"), "w").write(txt + "\n")
