"""Iteration 106: is the hemisphere split real or a survey effect? (PREREG_106.md)"""
import numpy as np, pandas as pd, os, gzip, pickle, warnings
from scipy.integrate import cumulative_trapezoid
from astropy.coordinates import SkyCoord; import astropy.units as u
warnings.filterwarnings("ignore")
here = os.path.dirname(os.path.abspath(__file__)); NAX = int(os.environ.get("NAX", 300))
zg = np.linspace(0, 2.5, 5001)
OMS = np.linspace(0.05, 0.8, 76)
DMF = [cumulative_trapezoid(1/np.sqrt(o*(1+zg)**3 + 1 - o), zg, initial=0) for o in OMS]
def uv(l, b): l, b = np.radians(l), np.radians(b); return np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)])
AX = uv(260, 12); rng = np.random.default_rng(5); RAX = [v/np.linalg.norm(v) for v in rng.normal(size=(NAX, 3))]
out = []
# ---------------- Check 1: Pantheon+ survey make-up and drop-one ----------------
PP = "/home/claude/pantheonplussh0es/datarelease/Pantheon+_Data/4_DISTANCES_AND_COVAR/"
sn = pd.read_csv(PP + "Pantheon+SH0ES.dat", sep=r"\s+")
raw = np.loadtxt(PP + "Pantheon+SH0ES_STAT+SYS.cov"); N = int(raw[0]); cov = raw[1:].reshape(N, N)
m = (sn["zHD"] > 0.05).values; sn = sn[m].reset_index(drop=True); cov = cov[np.ix_(m, m)]
g = SkyCoord(ra=sn.RA.values*u.deg, dec=sn.DEC.values*u.deg).galactic
lr, br = np.radians(g.l.deg), np.radians(g.b.deg); nh = np.c_[np.cos(br)*np.cos(lr), np.cos(br)*np.sin(lr), np.sin(br)]
z, zh, mb, sv = sn.zHD.values, sn.zHEL.values, sn.m_b_corr.values, sn.IDSURVEY.values
def om_fit(sel):
    Ci = np.linalg.inv(cov[np.ix_(sel, sel)]); cs = []
    for D in DMF:
        r = mb[sel] - 5*np.log10((1 + zh[sel])*np.interp(z[sel], zg, D)); B = (Ci @ r).sum(); cs.append(r @ Ci @ r - B**2/Ci.sum())
    cs = np.array(cs); i = np.argmin(cs); ok = OMS[cs <= cs[i] + 1]; return OMS[i], max((ok.max() - ok.min())/2, 0.005)
def split(keep, ax):
    t = (nh @ ax > 0) & keep; a = (nh @ ax <= 0) & keep
    (o1, e1), (o2, e2) = om_fit(t), om_fit(a); return o1 - o2, np.hypot(e1, e2), t.sum(), a.sum()
tow = nh @ AX > 0
out.append("CHECK 1 -- Pantheon+ (z > 0.05) survey make-up (IDSURVEY: count toward / away)")
ids, cnt = np.unique(sv, return_counts=True)
out.append("   " + ", ".join(f"{i}: {np.sum(tow & (sv == i))}/{np.sum(~tow & (sv == i))}" for i in ids[np.argsort(-cnt)]))
big = [i for i, c in zip(ids, cnt) if c > 0.05*len(sv)]
for drop in [None] + big:
    keep = np.ones(len(sv), bool) if drop is None else (sv != drop)
    d, e, nt, na = split(keep, AX)
    rd = np.array([split(keep, v)[0] for v in RAX]); p = np.mean(np.abs(rd) >= abs(d))
    out.append(f"   {'all surveys' if drop is None else f'drop survey {drop}':16s} N {nt}/{na}: matter-share difference (toward - away) {d:+.3f} +/- {e:.3f}; random-axis p = {p:.1%}")
# ---------------- Check 2: Union3 per-SN inputs, simplified standardisation ----------------
open(os.path.join(here, "iter106_checks.txt"), "w").write("\n".join(out) + "\n"); print("\n".join(out), flush=True)
def parse_coords(RA, DEC):
    ra, de = np.full(len(RA), np.nan), np.full(len(RA), np.nan)
    for i, (r, d) in enumerate(zip(RA, DEC)):
        try:
            r, d = str(r).strip(), str(d).strip()
            if ":" in r or " " in r:
                c = SkyCoord(r.replace(" ", ":"), d.replace(" ", ":"), unit=(u.hourangle, u.deg)); ra[i], de[i] = c.ra.deg, c.dec.deg
            else: ra[i], de[i] = float(r), float(d)
        except Exception: pass
    return ra, de
f = [x for x in os.listdir("/home/claude/union3_release") if x.endswith(".pickle")][0]
U = pickle.load(gzip.open("/home/claude/union3_release/" + f, "rb"), encoding="latin1")[0]
zc, mB, x1, c = (np.asarray(U[k], float) for k in ("z_CMB_list", "mB_list", "x1_list", "c_list"))
C3 = np.asarray(U["mBx1c_cov_list"], float); samp = np.asarray(U["sample_list"]); ra, de = parse_coords(U["RA"], U["Dec"])
sel = (zc > 0.05) & np.isfinite(ra) & np.isfinite(de) & (np.abs(c) < 0.3) & (np.abs(x1) < 3)
zc, mB, x1, c, C3, samp, ra, de = zc[sel], mB[sel], x1[sel], c[sel], C3[sel], samp[sel], ra[sel], de[sel]
g = SkyCoord(ra=ra*u.deg, dec=de*u.deg).galactic
lr, br = np.radians(g.l.deg), np.radians(g.b.deg); nu = np.c_[np.cos(br)*np.cos(lr), np.cos(br)*np.sin(lr), np.sin(br)]
smp = np.unique(samp); OH = (samp[:, None] == smp[None, :]).astype(float)
MU = [5*np.log10((1 + zc)*np.interp(zc, zg, D)) for D in DMF]
def chi(alpha, beta, s, k):
    var = C3[:, 0, 0] + alpha**2*C3[:, 1, 1] + beta**2*C3[:, 2, 2] + 2*alpha*C3[:, 0, 1] - 2*beta*C3[:, 0, 2] - 2*alpha*beta*C3[:, 1, 2] + 0.10**2
    w = 1/var[s]; y = (mB + alpha*x1 - beta*c - MU[k])[s]; X = OH[s][:, OH[s].sum(0) > 0]
    p = np.linalg.lstsq(X*np.sqrt(w)[:, None], y*np.sqrt(w), rcond=None)[0]; r = y - X @ p; return np.sum(w*r*r)
alls = np.ones(len(zc), bool)
best = min(((chi(a, b, alls, k), a, b) for a in np.linspace(0.10, 0.20, 11) for b in np.linspace(2.4, 3.6, 13) for k in range(0, 76, 5)))
_, AL, BE = best
def om_u(s):
    cs = np.array([chi(AL, BE, s, k) for k in range(76)]); i = np.argmin(cs); ok = OMS[cs <= cs[i] + 1]
    return OMS[i], max((ok.max() - ok.min())/2, 0.005)
def split_u(ax):
    t = nu @ ax > 0; (o1, e1), (o2, e2) = om_u(t), om_u(~t); return o1 - o2, np.hypot(e1, e2), t.sum(), (~t).sum(), o1, o2
d, e, nt, na, o1, o2 = split_u(AX)
rd = np.array([split_u(v)[0] for v in RAX]); p = np.mean(np.abs(rd) >= abs(d))
out += ["", f"CHECK 2 -- Union3 inputs (simplified: alpha {AL:.2f}, beta {BE:.2f}, per-sample offsets, 0.10 mag scatter), z > 0.05, N = {len(zc)}",
        f"   matter share toward {o1:.3f}, away {o2:.3f} (dark energy {1-o1:.3f} vs {1-o2:.3f}); difference {d:+.3f} +/- {e:.3f}; random-axis p = {p:.1%}  (N {nt}/{na})"]
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter106_checks.txt"), "w").write(txt + "\n")
print("Union3 coordinates parsed:", np.isfinite(parse_coords(U["RA"], U["Dec"])[0]).sum(), "of", len(U["RA"]))
