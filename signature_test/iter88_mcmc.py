import os, json, numpy as np
os.chdir("/home/claude/griduniverse/cosmology_fits")
exec(open("fit_law.py").read().split("Ndata = len(bao)")[0])
rng = np.random.default_rng(int.from_bytes(SNSET.encode(), "little") % 2**32)
def lnp(p):
    Om, h, wb, w0, wa = p
    if not (-2 < w0 < 0 and -3 < wa < 2): return -np.inf
    c = chi2("w0wa", p)
    return -0.5*c if c < 1e11 else -np.inf
x = np.array([0.31, 0.68, 0.0224, -0.8, -0.7]); lx = lnp(x)
step = np.array([0.008, 0.008, 0.00015, 0.06, 0.25])
chain = []
N = 24000
for i in range(N):
    y = x + step*rng.normal(size=5)*0.6
    ly = lnp(y)
    if np.log(rng.random()) < ly - lx: x, lx = y, ly
    if i > 4000 and i % 4 == 0: chain.append(x.copy())
chain = np.array(chain)
def zw(w0, wa):
    if wa == 0: return np.nan
    one_minus_a = -(1 + w0)/wa
    a = 1 - one_minus_a
    return 1/a - 1 if 0 < a <= 1 else np.nan
def zq(Om, h, w0, wa):
    Or = W_R/h**2; OL = 1 - Om - Or
    z = np.linspace(0, 3, 3001); zp = 1 + z; a = 1/zp
    rde = OL*zp**(3*(1 + w0 + wa))*np.exp(-3*wa*z/zp); w = w0 + wa*(1 - a)
    q = (Om*zp**3/2 + Or*zp**4 + rde*(1 + 3*w)/2)/(Om*zp**3 + Or*zp**4 + rde)
    s = np.where(np.diff(np.sign(q)) > 0)[0]
    return z[s[0]] if len(s) else np.nan
ZW = np.array([zw(p[3], p[4]) for p in chain]); ZQ = np.array([zq(p[0], p[1], p[3], p[4]) for p in chain])
ok = np.isfinite(ZW) & np.isfinite(ZQ) & (ZW < 3)
D = ZW[ok] - ZQ[ok]
res = dict(n=len(chain), frac_cross=float(ok.mean()), w0=[float(np.mean(chain[:,3])), float(np.std(chain[:,3]))],
           wa=[float(np.mean(chain[:,4])), float(np.std(chain[:,4]))], zw=[float(np.median(ZW[ok])), float(np.std(ZW[ok]))],
           zq=[float(np.median(ZQ[ok])), float(np.std(ZQ[ok]))], delta_median=float(np.median(D)), delta_mean=float(np.mean(D)),
           delta_std=float(np.std(D)), delta_16_84=[float(np.percentile(D, 16)), float(np.percentile(D, 84))],
           p_delta_gt0=float((D > 0).mean()))
print("RESULT", SNSET, json.dumps(res))
