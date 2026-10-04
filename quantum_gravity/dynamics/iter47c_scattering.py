"""
ITERATION 47c (follow-up, written after 47b showed a random grid scatters ~29% of a wavelength-10 wave in 16 time units and slows
the coherent wave by 2.4%; expectations committed BEFORE this run):
How do scattering and slowing scale with wavelength? This decides whether they matter for real gravitational waves.
 S1 expectation: point-disorder (Rayleigh-like) scattering -> loss per unit distance ~ k^4, i.e. slope ~ 4 in log-log.
 S2 if instead the slope is ~2 (gradient-noise scattering), the real-universe loss over 40 Mpc could approach ~1e-15 and must
    be checked against data; slope >= 3 makes it utterly negligible.
 S3 the coherent-speed deficit shrinks with wavelength at least as fast as k^2.
Method: same grid/operator as 47b; packets at wavelengths 7, 10, 14 (envelope width = wavelength/2); coherent energy and coherent
centroid measured as in 47b revision 2; compared with exact continuum evolution of the same launch.
"""
import numpy as np, scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.integrate import quad
rng = np.random.default_rng(470)
L, s, Rc = 40.0, 1.0, 3.5
n = int(L**3); X = rng.uniform(0, L, (n, 3)); tree = cKDTree(X, boxsize=L)
Dm = tree.sparse_distance_matrix(tree, Rc, output_type="coo_matrix")
W = Dm.copy(); W.data = np.exp(-Dm.data**2/(2*s**2)); W = W.tocsr(); W.setdiag(0); W.eliminate_zeros()
norm = quad(lambda r: 4*np.pi*r**4/3*np.exp(-r**2/(2*s**2)), 0, Rc)[0]/2
Lop = (W - sp.diags(np.asarray(W.sum(1)).ravel()))/norm
w2 = np.vectorize(lambda k: quad(lambda r: 4*np.pi*r**2*np.exp(-r**2/2)*(1 - np.sinc(k*r/np.pi)), 0, Rc)[0]/norm)
xs = X[:, 0]; bins = np.linspace(0, L, 161); bc = 0.5*(bins[1:] + bins[:-1]); cnt = np.histogram(xs, bins)[0]
xg = np.linspace(0, L, 2048, endpoint=False); kg = 2*np.pi*np.fft.fftfreq(2048, L/2048); omg = np.sqrt(w2(np.abs(kg)))
def run(lam, T=16.0, dt=0.05, x0=10.0):
    k0 = 2*np.pi/lam; sig = lam/2; wk = np.sqrt(w2(k0))
    def launch(x):
        d = (x - x0 + L/2) % L - L/2; g = np.exp(-d**2/(2*sig**2))
        return g*np.cos(k0*x), -(wk/k0)*((-d/sig**2)*g*np.cos(k0*x) - k0*g*np.sin(k0*x))
    h, v = launch(xs)
    def coh(h, v):
        hb = np.histogram(xs, bins, weights=h)[0]/cnt; vb = np.histogram(xs, bins, weights=v)[0]/cnt
        return hb**2 + vb**2/wk**2
    E0 = coh(h, v).sum(); cur = x0; ts, cs, fr = [], [], []
    v = v + 0.5*dt*(Lop @ h)
    for i in range(int(T/dt) + 1):
        if i % 20 == 0:
            vt = v - 0.5*dt*(Lop @ h); e = coh(h, vt); d = (bc - cur + L/2) % L - L/2; m = np.abs(d) < 2*lam
            cur = cur + np.sum(e[m]*d[m])/np.sum(e[m]); ts.append(i*dt); cs.append(cur); fr.append(e.sum()/E0)
        h = h + dt*v; v = v + dt*(Lop @ h)
    vgrid = np.polyfit(ts, cs, 1)[0]
    # continuum, same launch
    F, V = launch(xg); Fh, Vh = np.fft.fft(F), np.fft.fft(V); oms = np.where(omg > 0, omg, 1); cc = []; cur = x0
    for t in ts:
        hh = np.real(np.fft.ifft(Fh*np.cos(omg*t) + Vh*np.sin(omg*t)/oms)); vv = np.real(np.fft.ifft(-Fh*omg*np.sin(omg*t) + Vh*np.cos(omg*t)))
        e = hh**2 + vv**2/wk**2; d = (xg - cur + L/2) % L - L/2; m = np.abs(d) < 2*lam; cur = cur + np.sum(e[m]*d[m])/np.sum(e[m]); cc.append(cur)
    vcont = np.polyfit(ts, cc, 1)[0]
    rate = -np.polyfit(np.array(ts)*vcont, np.log(fr), 1)[0]          # coherent energy loss per unit distance
    return k0, vgrid/vcont, rate
out = ["ITERATION 47c: scattering of waves by the random grid vs wavelength (expectations committed first)", ""]
res = []
for lam in (7.0, 10.0, 14.0):
    k0, ratio, rate = run(lam); res.append((k0, ratio, rate))
    out.append(f"wavelength {lam:4.1f}: coherent speed / continuum = {ratio:.4f}; coherent energy lost per unit distance = {rate:.4f}")
k = np.array([r[0] for r in res]); rt = np.array([r[2] for r in res]); dv = 1 - np.array([r[1] for r in res])
p = np.polyfit(np.log(k), np.log(rt), 1)[0]; out.append(f"loss-rate slope d ln(rate)/d ln k = {p:.2f}")
if np.all(dv > 0): out.append(f"speed-deficit slope = {np.polyfit(np.log(k), np.log(dv), 1)[0]:.2f}")
txt = "\n".join(out); print(txt); open("iter47c_scattering.txt", "w").write(txt + "\n")
