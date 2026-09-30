"""
Growth-rate test (f sigma8) of the dark-energy law against LCDM and w0wa.
Predictions: each model's best fit to Planck CMB + DESI DR2 BAO + supernovae (the Cobaya/CAMB minima), run through CAMB.
Data: SDSS DR16 consensus redshift-space-distortion growth measurements (Alam et al. 2021), full covariance,
      marginalised over the distance parameters:  BOSS z = 0.38, 0.51 ; eBOSS LRG z = 0.698 ; eBOSS ELG z = 0.845 (likelihood grid) ;
      eBOSS QSO z = 1.48.   This is an out-of-sample test: no growth data went into the fits.
"""
import numpy as np, camb, glob, json
from scipy.integrate import solve_ivp
D = "/home/claude/cobayasampler/bao_data/"
# ---------- data ----------
cl = np.loadtxt(D + "sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8_covtot.txt"); cq = np.loadtxt(D + "sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8_covtot.txt")
vl = np.array([float(l.split()[1]) for l in open(D + "sdss_DR16_BAOplus_LRG_FSBAO_DMDHfs8.dat") if l.strip()])
vq = np.array([float(l.split()[1]) for l in open(D + "sdss_DR16_BAOplus_QSO_FSBAO_DMDHfs8.dat") if l.strip()])
iL = [2, 5, 8]
z_g = np.array([0.38, 0.51, 0.698, 1.48]); fs_g = np.r_[vl[iL], vq[2]]
C = np.zeros((4, 4)); C[:3, :3] = cl[np.ix_(iL, iL)]; C[3, 3] = cq[2, 2]; Ci = np.linalg.inv(C)
g = np.loadtxt(D + "sdss_DR16_ELG_FSBAO_DMDHfs8gridlikelihood.txt")
fsv = np.unique(g[:, 2]); prob = np.array([g[g[:, 2] == f, 3].sum() for f in fsv]); prob /= prob.max()
def chi2_elg(fs): return -2 * np.log(max(np.interp(fs, fsv, prob), 1e-300))
cen = fsv[np.argmax(prob)]; m1 = np.sum(fsv * prob) / prob.sum(); s1 = np.sqrt(np.sum((fsv - m1)**2 * prob) / prob.sum())
print("data  z:", list(z_g) + [0.845], " f sigma8:", np.round(fs_g, 3).tolist() + [round(m1, 3)], " sd:", np.round(np.sqrt(np.diag(C)), 3).tolist() + [round(s1, 3)])

# ---------- models ----------
def law_w(H0, ombh2, omch2, omnuh2, beta=0.5):
    h = H0 / 100; om = (ombh2 + omch2 + omnuh2) / h**2; orad = 4.18e-5 / h**2; ol = 1 - om - orad
    lg = np.linspace(0.0, np.log(1e-5), 600)
    rhs = lambda lna, y: [beta * (om*np.exp(-3*lna)/2 + orad*np.exp(-4*lna) - np.exp(y[0])) / (om*np.exp(-3*lna) + orad*np.exp(-4*lna) + np.exp(y[0]))]
    s = solve_ivp(rhs, [0, lg[-1]], [np.log(ol)], t_eval=lg, rtol=1e-8, atol=1e-10)
    a = np.exp(s.t); rde = np.exp(s.y[0]); rm, rr = om*a**-3, orad*a**-4; q = (rm/2 + rr - rde)/(rm + rr + rde)
    return a[::-1], (-1 - beta*q/3)[::-1]

def best_point(pattern):
    best = None
    for f in glob.glob(pattern):
        hdr = open(f).readline().lstrip("#").split(); v = np.loadtxt(f); d = dict(zip(hdr, v))
        if best is None or d["chi2"] < best["chi2"]: best = d; best["file"] = f.split("/")[-1]
    return best

zs = np.array([0.0, 0.38, 0.51, 0.698, 0.845, 1.48])
def predict(model, p):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], tau=p["tau"], mnu=0.06, nnu=3.044, num_massive_neutrinos=1)
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    if model == "w0wa": pars.set_dark_energy(w=p["w"], wa=p["wa"], dark_energy_model="ppf")
    if model == "LAW":
        a, w = law_w(p["H0"], p["ombh2"], p["omch2"], pars.omnuh2); pars.set_dark_energy_w_a(a, w, dark_energy_model="ppf")
    pars.set_matter_power(redshifts=zs[::-1].tolist(), kmax=2.0); pars.WantTransfer = True
    r = camb.get_results(pars)
    fs8 = r.get_fsigma8()[::-1]; s8 = r.get_sigma8()[::-1]
    Om = (p["ombh2"] + p["omch2"] + pars.omnuh2) / (p["H0"]/100)**2
    return fs8, s8, Om

out = {}
for sn, tag in [("DES-Dovekie", "desdovekie"), ("Pantheon+", "pantheonplus")]:
    print(f"\n=== fits using {sn} supernovae ===")
    rows = {}
    for model, pat in [("LCDM", f"LCDM_sn.{tag}*.minimum.txt"), ("LAW", f"LAW_0.5_sn.{tag}*.minimum.txt"), ("w0wa", f"w0wa_sn.{tag}*.minimum.txt")]:
        p = best_point("/home/claude/fullfit/out/" + pat)
        fs8, s8, Om = predict(model, p)
        pred = fs8[[1, 2, 3, 5]]; d = pred - fs_g; c2 = d @ Ci @ d + chi2_elg(fs8[4])
        # growth amplitude preferred by the data relative to this model: fs_data = A * fs_model
        A = (pred @ Ci @ fs_g) / (pred @ Ci @ pred); Aerr = 1 / np.sqrt(pred @ Ci @ pred)
        S8 = s8[0] * np.sqrt(Om / 0.3)
        rows[model] = dict(chi2_growth=float(c2), fs8=fs8.tolist(), sigma8=float(s8[0]), S8=float(S8), A=float(A), Aerr=float(Aerr), file=p["file"], chi2_main=float(p["chi2"]))
        print(f"{model:5s} ({p['file']}): f sigma8 at z = 0.38/0.51/0.70/0.85/1.48 = {np.round(fs8[1:], 3).tolist()} | sigma8 = {s8[0]:.3f}, S8 = {S8:.3f} | "
              f"growth chi2 = {c2:.2f} (5 points) | data/model growth amplitude = {A:.3f} ± {Aerr:.3f}")
    for m in ["LAW", "w0wa"]:
        print(f"   {m} - LCDM: growth delta chi2 = {rows[m]['chi2_growth'] - rows['LCDM']['chi2_growth']:+.2f};  combined with CMB+BAO+SN: {rows[m]['chi2_growth'] + rows[m]['chi2_main'] - rows['LCDM']['chi2_growth'] - rows['LCDM']['chi2_main']:+.2f}")
    out[sn] = rows
json.dump(dict(data=dict(z=z_g.tolist() + [0.845], fs8=fs_g.tolist() + [m1], sd=np.sqrt(np.diag(C)).tolist() + [s1]), fits=out),
          open("/home/claude/fullfit/growth_results.json", "w"), indent=1)
