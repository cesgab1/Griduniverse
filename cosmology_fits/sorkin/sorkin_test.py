"""Sorkin's everpresent Λ (causal-set grid): Λ = α S(V)/V, S a random walk in the 4-volume V of the causal past (Planck length cancels).
Each realisation is evolved from the radiation era to today (T_CMB fixed; ω_b, ω_c fixed from the CMB); H0 comes out.
Scored with DESI DR2 BAO + DES-Dovekie SN (cobaya likelihoods) + the CMB acoustic angle θ*; compared with ΛCDM and the β = ½ law (H0 fitted)."""
import sys, numpy as np, json
from scipy.optimize import minimize_scalar
sys.path.insert(0, "/home/claude/fullfit/sorkin")
from cobaya.model import get_model
import bgtheory
Mpc=3.0857e22; c=2.998e8; G=6.674e-11; Ckm=299792.458
wb, wc, wg = 0.02237, 0.1200, 2.47e-5; wr = wg*(1+0.2271*3.044)
rc = 3*(100e3/Mpc)**2/(8*np.pi*G)
rho_m0, rho_r0 = (wb+wc)*rc, wr*rc
zstar, zdrag = 1089.9, 1059.9
CNT = 0
model = get_model({"packages_path": "/home/claude/cobaya_packages", "likelihood": {"bao.desi_dr2": None, "sn.desdovekie": None},
                   "theory": {"bgtheory.BG": None}, "params": {"dummy": {"prior": {"min": 0, "max": 1e12}}, "rdrag": {"derived": True}}}, stop_at_error=True)
zgrid = np.concatenate([np.linspace(0, 3, 600), np.geomspace(3.01, 1500, 400)])
def horizon(lna, H, zend):
    a = np.exp(lna); R = 3*wb/(4*wg)*a; cs = c/np.sqrt(3*(1+R)); m = a <= 1/(1+zend)
    return np.trapezoid((cs/(a*H))[m], lna[m])/Mpc
def raw(lna, H):
    a = np.exp(lna); zz = 1/a-1; o = np.argsort(zz); Hz = np.interp(zgrid, zz[o], H[o]*Mpc/1e3)
    return Hz, horizon(lna, H, zdrag), horizon(lna, H, zstar)
LNA = np.linspace(np.log(1e-8), 0, 3000)
def lcdm_H(h, beta=0.0):
    a = np.exp(LNA); om, orad = (wb+wc)/h**2, wr/h**2; ol = 1-om-orad
    if beta == 0: rde = ol*np.ones_like(a)
    else:
        from scipy.integrate import solve_ivp
        f = lambda x, y: [beta*(om*np.exp(-3*x)/2 + orad*np.exp(-4*x) - np.exp(y[0]))/(om*np.exp(-3*x)+orad*np.exp(-4*x)+np.exp(y[0]))]
        s = solve_ivp(f, [0, LNA[0]], [np.log(ol)], t_eval=LNA[::-1], rtol=1e-8, atol=1e-10); rde = np.exp(s.y[0][::-1])
    return h*100e3/Mpc*np.sqrt(om*a**-3 + orad*a**-4 + rde)
_, RD0, RS0 = raw(LNA, lcdm_H(0.674))           # reference normalisation: CAMB gives r_d 147.09, r* 144.39 for this cosmology
def score(lna, H):
    Hz, rd, rs = raw(lna, H); rd *= 147.09/RD0; rs *= 144.39/RS0
    bgtheory.BG.state_bg = {"z": zgrid, "H": Hz, "rdrag": rd}
    global CNT; CNT += 1
    ll = model.loglikes({"dummy": float(CNT)}, as_dict=True)[0]
    chi = np.concatenate([[0], np.cumsum(0.5*(Ckm/Hz[1:]+Ckm/Hz[:-1])*np.diff(zgrid))])
    th = 100*rs/np.interp(zstar, zgrid, chi)
    L = model.likelihood["bao.desi_dr2"]
    thv = np.array([L.theory_fun(z, o) for z, o in zip(L.data["z"], L.data["observable"])]).ravel(); rr = L.data["value"].values - thv
    d = {"bao": float(rr @ L.invcov @ rr), "sn": -2*ll["sn.desdovekie"], "theta": ((th-1.04109)/0.0003)**2, "H0": Hz[0]}
    d["tot"] = d["bao"] + d["sn"] + d["theta"]; return d
def best(beta):
    r = minimize_scalar(lambda h: score(LNA, lcdm_H(h, beta))["tot"], bounds=(0.6, 0.75), method="bounded")
    return score(LNA, lcdm_H(r.x, beta))
def evolve(alpha, seed, n=2200):
    rng = np.random.default_rng(seed)
    lna = np.linspace(np.log(1e-8), 0, n); a = np.exp(lna); dl = lna[1]-lna[0]
    H = np.zeros(n); t = np.zeros(n); eta = np.zeros(n); S = 0.0; Vp = 0.0; Lam = 0.0; dt = np.zeros(n)
    for i in range(n):
        rho = rho_m0*a[i]**-3 + rho_r0*a[i]**-4
        H2 = 8*np.pi*G/3*rho + Lam*c**2/3
        if H2 <= 0: return None, None
        H[i] = np.sqrt(H2)
        if i == 0: t[0] = 1/(2*H[0]); eta[0] = 2*t[0]/a[0]; dt[0] = t[0]
        else: dt[i] = dl/H[i]; t[i] = t[i-1]+dt[i]; eta[i] = eta[i-1]+dl/(a[i]*H[i])
        V = np.sum(a[:i+1]**3*4*np.pi/3*(c*(eta[i]-eta[:i+1]))**3*c*dt[:i+1])
        S += np.sqrt(max(V-Vp, 0.0))*rng.standard_normal(); Vp = V
        Lam = alpha*S/V if (V > 0 and i >= 150) else 0.0
    return lna, H
if __name__ == "__main__":
    out = {"lcdm": best(0.0), "beta_half": best(0.5)}
    print("ΛCDM", out["lcdm"]); print("β=½ ", out["beta_half"])
    for alpha in [float(x) for x in sys.argv[1].split(",")]:
        res = []
        for seed in range(int(sys.argv[2])):
            lna, H = evolve(alpha, seed)
            if lna is None: res.append(None); continue
            s = score(lna, H); s["OL0"] = 1 - (8*np.pi*G/3*(rho_m0+rho_r0))/H[-1]**2; res.append(s)
        ok = [r for r in res if r is not None]
        tot = np.array([r["tot"] for r in ok])
        if len(ok) == 0:
            print(f"α={alpha}: 0/{len(res)} survive (negative Λ -> universe recollapses in every realisation)"); out[f"alpha_{alpha}"] = res; continue
        print(f"α={alpha}: {len(ok)}/{len(res)} survive; Δχ² vs ΛCDM: best {tot.min()-out['lcdm']['tot']:+.1f}, "
              f"median {np.median(tot)-out['lcdm']['tot']:+.1f}; share beating ΛCDM {np.mean(tot<out['lcdm']['tot']):.1%}, beating β=½ {np.mean(tot<out['beta_half']['tot']):.1%}")
        out[f"alpha_{alpha}"] = res
    json.dump(out, open(f"sorkin_{sys.argv[1]}.json", "w"), default=float)
