"""
Round 2 of public-data checks (none of these data went into any fit):
 1. Cosmic chronometers: 32 direct H(z) measurements (Moresco et al. compilation; diagonal errors, systematics covariance not included)
 2. Local Hubble constant: SH0ES (Riess et al. 2022) 73.04 +- 1.04 ; CCHP TRGB (Freedman et al. 2025) 70.39 +- 1.94
 3. Age of the Universe from the oldest globular clusters (Valcin et al. 2025): 13.57 +- 0.27 Gyr
 4. Black-hole shadows (EHT): Sgr A* ring 51.8 +- 2.3 uas; M87* ring 42 +- 3 uas  vs the framework's photon-orbit shadow 2 sqrt(27) GM/(c^2 D)
 5. Speed of gravity (GW170817 + GRB 170817A): -3e-15 < (c_g - c)/c < 7e-16, vs the grid's dispersion at LIGO frequencies
Predictions: each model's CMB + BAO + SN full-likelihood best fit (Cobaya/CAMB minima) run through CAMB.
"""
import numpy as np, json, glob, camb
from scipy.integrate import solve_ivp
src = open("/home/claude/fullfit/growth_test.py").read().split("out = {}")[0]
src = src.split("# ---------- models ----------")[1]
exec(src)
cc = np.array([[0.070,69,19.6],[0.090,69,12],[0.120,68.6,26.2],[0.170,83,8],[0.1791,75,4],[0.1993,75,5],[0.200,72.9,29.6],[0.270,77,14],
 [0.280,88.8,36.6],[0.3519,83,14],[0.3802,83,13.5],[0.400,95,17],[0.4004,77,10.2],[0.4247,87.1,11.2],[0.445,92.8,12.9],[0.470,89,50],
 [0.4783,80.9,9],[0.480,97,62],[0.5929,104,13],[0.6797,92,8],[0.750,98.8,33.6],[0.7812,105,12],[0.8754,125,17],[0.880,90,40],[0.900,117,23],
 [1.037,154,20],[1.300,168,17],[1.363,160,33.6],[1.430,177,18],[1.530,140,14],[1.750,202,40],[1.965,186.5,50.4]])
def background(model, p):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], tau=p["tau"], mnu=0.06, nnu=3.044, num_massive_neutrinos=1)
    pars.InitPower.set_params(As=p["As"], ns=p["ns"])
    if model == "w0wa": pars.set_dark_energy(w=p["w"], wa=p["wa"], dark_energy_model="ppf")
    if model == "LAW":
        a, w = law_w(p["H0"], p["ombh2"], p["omch2"], pars.omnuh2); pars.set_dark_energy_w_a(a, w, dark_energy_model="ppf")
    return camb.get_background(pars)
out = {}
for sn, tag in [("DES-Dovekie", "desdovekie"), ("Pantheon+", "pantheonplus")]:
    print(f"\n=== best fits using {sn} supernovae ===")
    for model, pat in [("LCDM", f"LCDM_sn.{tag}*.minimum.txt"), ("LAW", f"LAW_0.5_sn.{tag}*.minimum.txt"), ("w0wa", f"w0wa_sn.{tag}*.minimum.txt")]:
        p = best_point("/home/claude/fullfit/out/" + pat); bg = background(model, p)
        Hz = bg.hubble_parameter(cc[:, 0]); c2cc = np.sum(((Hz - cc[:, 1]) / cc[:, 2])**2)
        H0 = p["H0"]; age = bg.physical_time(0)
        t_sh = (H0 - 73.04) / np.hypot(1.04, 0.5); t_cc = (H0 - 70.39) / np.hypot(1.94, 0.5)
        t_age = (age - 13.57) / np.hypot(0.27, 0.1)
        out[f"{sn}_{model}"] = dict(chi2_cc=float(c2cc), H0=float(H0), age=float(age), t_shoes=float(t_sh), t_cchp=float(t_cc), t_age=float(t_age))
        print(f"{model:5s}: chronometers chi2 = {c2cc:6.2f} / 32 | H0 = {H0:.2f} (vs SH0ES {t_sh:+.1f} sigma, vs CCHP {t_cc:+.1f} sigma) | "
              f"age = {age:.2f} Gyr (vs oldest stars {t_age:+.1f} sigma)")
G, c, Msun, kpc, Mpc = 6.6743e-11, 2.99792458e8, 1.98847e30, 3.0857e19, 3.0857e22
uas = np.pi / 180 / 3600 * 1e-6
print("\n=== Black-hole shadows (framework = GR outside the horizon; the capped core shifts it by ~(l_P/r)^2 ~ 1e-80) ===")
for name, M, D, ring, er in [("Sgr A* (GRAVITY mass 4.297e6, 8.277 kpc)", 4.297e6, 8.277 * kpc, 51.8, 2.3),
                             ("M87* (stellar-dynamics mass 6.2e9 at 16.8 Mpc)", 6.2e9, 16.8 * Mpc, 42.0, 3.0)]:
    th = 2 * np.sqrt(27) * G * M * Msun / (c**2 * D) / uas
    print(f"   {name}: predicted shadow {th:.1f} uas, EHT ring {ring} ± {er} uas, ring/shadow = {ring/th:.2f} (emission rings sit ~0-10% outside the shadow)")
f, ld = 100.0, 0.603 * 1.616255e-35
kl = 2 * np.pi * f / c * ld
print(f"\n=== Speed of gravity: grid dispersion at 100 Hz gives (c_g - c)/c ~ -0.0375 (k l_d)^2 = {-0.0375 * kl**2:.1e}; GW170817 allows -3e-15 to +7e-16 -> passes")
json.dump(out, open("/home/claude/fullfit/data_round2.json", "w"), indent=1)
