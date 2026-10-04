"""
ITERATION 21: size of dark energy from random SPIN flips of the grid's cells (Coalesce: spin as the coin), combined with
the grid's memory (kappa = 3 from Hubble friction).
Model: every Planck 4-cell in the remembered past region (memory 1/(kappa H) back in time, causal radius c/(kappa H))
contributes +1 or -1 at random. Leftover = sqrt(N) of N nudges -> tension per unit 4-volume ~ 1/sqrt(N) (Planck units).
   N(t) = (c/(kappa H))^3 (1/(kappa H)) / (l_P^3 t_P)  ->  rho_DE ~ alpha * rho_P * kappa^2 (l_P H)^2 * xi(t)
   with xi a unit random variable that is renewed every memory time (an OU process in time, rate kappa H), alpha = O(1).
Equivalent: Omega_DE(t) = c_alpha * xi(t), i.e. the SIZE always tracks the critical density (that is why it is the right
order today). This is Sorkin's everpresent-Lambda scaling with the toy's memory window. Checks:
 A. magnitude today for alpha = 1
 B. history: Omega_DE must be tiny early (CMB: early dark energy <~ 0.02 near recombination; BBN <~ 0.1) and ~0.69 today.
    Monte Carlo of xi over ln a with correlation 1/kappa e-fold; scan c_alpha; probability that a history passes both.
"""
import numpy as np
rng = np.random.default_rng(21)
lP, tP = 1.616e-35, 5.39e-44; H0 = 67.4e3/3.0857e22; c = 2.998e8; kappa = 3.0
rhoP = 5.16e96; rho_crit = 3*H0**2/(8*np.pi*6.674e-11)
N = (c/(kappa*H0))**3*(1/(kappa*H0))/(lP**3*tP)
rho1 = rhoP/np.sqrt(N)
out = ["ITERATION 21: dark-energy size from random spin flips of grid cells, with the toy's memory", "",
       f"A. nudges remembered now: N = {N:.1e}; leftover density (alpha = 1) = {rho1:.1e} kg/m^3 = {rho1/rho_crit:.0f} x critical",
       f"   (observed dark energy = 0.69 x critical): right order of magnitude - the size problem's 120 orders become a factor ~{rho1/rho_crit/0.69:.0f}",
       "", "B. but the size then TRACKS the critical density at every epoch: Omega_DE(t) = c_alpha * xi(t), xi random (memory kappa = 3)"]
lna = np.linspace(np.log(1/1e9), 0, 6000); d = lna[1] - lna[0]                 # from BBN (z ~ 1e9) to today
nrun = 20000; xi = np.zeros((nrun, len(lna))); xi[:, 0] = rng.standard_normal(nrun)
rho_ = np.exp(-kappa*d); s = np.sqrt(1 - rho_**2)
for i in range(1, len(lna)): xi[:, i] = rho_*xi[:, i - 1] + s*rng.standard_normal(nrun)
z = np.exp(-lna) - 1; rec = (z > 800) & (z < 1500); bbn = z > 1e8
for ca in (0.005, 0.01, 0.02, 0.05, 0.1, 0.3, 1.0):
    Om = ca*xi; ok_early = (np.abs(Om[:, rec]).max(1) < 0.02) & (np.abs(Om[:, bbn]).max(1) < 0.1)
    ok_today = (Om[:, -1] > 0.62) & (Om[:, -1] < 0.76)
    out.append(f"   c_alpha = {ca:5.3f}: early limits pass {ok_early.mean():6.1%}; today in 0.62-0.76: {ok_today.mean():8.3%}; BOTH: {(ok_early & ok_today).mean():.4%}")
out += ["", "C. recent history even in the passing cases: Omega_DE jumps every ~1/3 e-fold (memory), so w(z) is noisy, not the smooth",
        "   single crossing the data prefer (Claim 1 / toy fits); and half of all draws are NEGATIVE.",
        "", "Reading: the random spin flips give the right SIZE today, but a size that is recomputed from the remembered region at",
        "every epoch makes dark energy a constant FRACTION of the universe (plus noise). The CMB and nucleosynthesis need it",
        "negligible early, so today's 0.69 would be an extremely unlikely draw. The size and the shape are not independent:",
        "the size must be accumulated or frozen, not recomputed from the horizon each moment."]
txt = "\n".join(out); print(txt); open("iter21_spin_nudges.txt", "w").write(txt + "\n")
