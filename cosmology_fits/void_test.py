"""
Fast-motion (slack grid / MOND-type) picture -> we sit in a large emptied region (void); matter streaming out of it
makes nearby expansion look faster. Test against Pantheon+SH0ES (Cepheid-calibrated, full covariance, all z > 0.0233).
Void outflow: peculiar velocity v(d) = eps * H_out * d / sqrt(1 + (d/R)^6)   (uniform excess expansion eps inside radius R,
falling as 1/d^2 outside, like an underdense region), observer near the centre.
Observed redshift (1+z) = (1+z_cos(d))(1+v/c);  d_L = d_L,cos(d) (1+v/c)^2.
H_out (the universe beyond the void) is what the early-universe ruler measures: Planck + DESI give 67.9 +- 0.5.
Redshifts: zCMB (no flow corrections, so a real outflow is not partly removed); zHD also shown.
"""
import numpy as np, pandas as pd, json, sys
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
D = "/home/claude/cobaya_packages/data/sn_data/PantheonPlus/"
d = pd.read_csv(D + "Pantheon+SH0ES.dat", sep=r"\s+"); n = len(d)
C = np.loadtxt(D + "Pantheon+SH0ES_STAT+SYS.cov", skiprows=1).reshape(n, n)
c = 299792.458; Om = float(sys.argv[1]) if len(sys.argv) > 1 else 0.315
cal = d.IS_CALIBRATOR.values == 1
zg = np.linspace(0, 2.5, 20001); Ez = np.sqrt(Om*(1+zg)**3 + 1 - Om); chi_g = cumulative_trapezoid(1/Ez, zg, initial=0)  # in c/H units
out = {}
for zcol in ("zCMB", "zHD"):
    z = d[zcol].values; zh = d.zHEL.values
    use = cal | (z > 0.0233); idx = np.where(use)[0]; Ci = np.linalg.inv(C[np.ix_(idx, idx)]); one = np.ones(len(idx))
    mobs = d.m_b_corr.values[idx]; iscal = cal[idx]; zz = z[idx]; zhh = zh[idx]; ceph = d.CEPH_DIST.values[idx]
    def mu_model(Hout, eps, R):
        dcom = chi_g * c / Hout                                   # comoving distance grid (Mpc)
        v = eps * Hout * dcom / np.sqrt(1 + (dcom / R)**6) if eps else 0*dcom
        zobs = (1 + zg) * (1 + v/c) - 1
        if np.any(np.diff(zobs) <= 0): return None
        dcos = np.interp(zz, zobs, dcom); zc = np.interp(zz, zobs, zg); vv = np.interp(zz, zobs, v)
        dL = dcos * (1 + zc) * (1 + vv/c)**2 * (1 + zhh)/(1 + zz)
        return 5*np.log10(dL) + 25
    def chi2(Hout, eps, R, prior=True):
        mu = mu_model(Hout, eps, R)
        if mu is None: return 1e9
        r = mobs - np.where(iscal, ceph, mu)
        M = (one @ Ci @ r)/(one @ Ci @ one); r = r - M          # brightness scale M marginalised (best fit)
        return r @ Ci @ r + (((Hout - 67.9)/0.5)**2 if prior else 0)
    res = {}
    # A: no void, H0 free (supernovae + Cepheids only)
    a = minimize(lambda p: chi2(p[0], 0, 1, False), [73], method="Nelder-Mead"); res["A_noVoid_free"] = (a.fun, a.x[0])
    # B: no void, H0 = early-universe value
    res["B_noVoid_ruler"] = (chi2(67.9, 0, 1, False) , 67.9)
    # C: void, H_out tied to the ruler, eps and R free (scan R)
    best = None; scan = []
    for R in (100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400):
        b = minimize(lambda p: chi2(p[0], p[1], R), [67.9, 0.06], method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3))
        scan.append((R, b.fun, b.x[0], b.x[1]))
        if best is None or b.fun < best[1]: best = (R, b.fun, b.x[0], b.x[1])
    res["C_void_scan"] = scan; res["C_best"] = best
    # D: the Mazurenko et al.-type void: about half density out to 300 Mpc, radius ~500 Mpc; linear-theory outflow eps ~ f*|delta|/3 (MOND faster growth): take eps from data at R=500
    b = minimize(lambda p: chi2(p[0], p[1], 500), [67.9, 0.07], method="Nelder-Mead"); res["D_R500"] = (b.fun, b.x[0], b.x[1])
    out[zcol] = res
    print(f"\n=== redshifts {zcol}, Om = {Om}, {len(idx)} light curves ===")
    print(f"A  no void, H0 free (supernovae only):          chi2 = {res['A_noVoid_free'][0]:.1f}   H0 = {res['A_noVoid_free'][1]:.2f}")
    print(f"B  no void, H0 = 67.9 (early-universe ruler):   chi2 = {res['B_noVoid_ruler'][0]:.1f}   (worse by {res['B_noVoid_ruler'][0]-res['A_noVoid_free'][0]:.1f}: the Hubble tension)")
    print(f"C  void, H_out = ruler (prior 67.9 +- 0.5), best over size:")
    for R, f, H, e in scan: print(f"     void radius {R:5d} Mpc: chi2 = {f:.1f}  H_out = {H:.2f}  local excess eps = {e*100:+.1f}%")
    print(f"   best: R = {best[0]} Mpc, chi2 = {best[1]:.1f}  -> tension relieved by {res['B_noVoid_ruler'][0]-best[1]:.1f} of {res['B_noVoid_ruler'][0]-res['A_noVoid_free'][0]:.1f}")
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/void_test.json", "w"), indent=1, default=float)
