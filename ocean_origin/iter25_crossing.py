"""
ITERATION 25 (pre-registration part): can 'crossing between sheets and gaps' explain the dark/ordinary ratio 5.36?
Picture (Coalesce question, Oct 2026): early on, while hot and dense, things could cross between the sheets (ordinary matter)
and the gaps (the Ocean). A conserved leftover charge (the matter-over-antimatter excess, 'baryon number') was SHARED. When
crossing shut at temperature T_f, the split froze. BORROWED framework: asymmetric dark matter (Kaplan, Luty, Zurek 2009).
What it predicts:
  number ratio n_Ocean/n_baryon = (charge capacity of the gap content)/(charge capacity of the sheet content) at T_f
  -> Ocean particle mass m = 5.36 m_p * n_B/n_Ocean           (the 5.36 becomes a MASS prediction, not an explanation)
  -> crossing rate = expansion rate at T_f fixes the coupling G (contact interaction: Gamma ~ G^2 T^5 = H(T_f))
  -> today's scattering off nuclei sigma ~ G^2 mu^2 / pi        (underground detectors)
  -> no annihilation signal today (no anti-Ocean left): consistent with Fermi/AMS/IceCube seeing nothing
Free choices (counted): (1) gap content = one Dirac fermion carrying charge 1 (simplest; scalar gives the same capacity);
  charge 1/3 shown as alternative; (2) contact interaction (heavy go-between); (3) T_f scanned 0.2-100 GeV (QCD to electroweak);
  electric-charge neutrality and sphaleron constraints ignored (O(1) effects, standard ADM papers).
Order-of-magnitude for sigma (factor ~10).
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
mp = 0.93827; R_obs = 0.1200/0.02237     # Planck omega_c / omega_b = 5.36

def cap(m, T, g, q, boson=False):
    """charge capacity: Delta n = q^2 g mu T^2/6 * c(m/T) in units where massless fermion gives c = 1 (boson 2 at m=0)"""
    x = m/T
    s = -1 if boson else 1
    f = lambda y: y**2*np.exp(-np.sqrt(y*y + x*x))/(1 + s*np.exp(-np.sqrt(y*y + x*x)))**2
    c = quad(f, 0, 60 + x)[0]/(np.pi**2/6)   # massless fermion integral = pi^2/6
    return g*q*q*c

quarks = [("u", 0.0022), ("d", 0.0047), ("s", 0.095), ("c", 1.27), ("b", 4.18), ("t", 173.0)]
def sheet_cap(T):
    return sum(cap(m, T, 6, 1/3) for _, m in quarks)   # 3 colours x 2 spins, baryon charge 1/3

def gstar(T):   # rough relativistic dof for H(T)
    return 10.75 if T < 0.15 else (61.75 if T < 1.0 else (75.75 if T < 4 else (86.25 if T < 80 else 106.75)))
MPl = 1.2209e19
out = ["ITERATION 25: sheet <-> gap crossing (asymmetric dark matter reading), predictions before checking detector limits", "",
       f"observed omega_c/omega_b = {R_obs:.3f}", "",
       "   T_f [GeV]   gap content          Ocean mass m [GeV]   coupling G [GeV^-2]   sigma_nucleon [cm^2]"]
rows = []
for qchi, label in ((1.0, "Dirac fermion, q=1"), (1/3, "Dirac fermion, q=1/3")):
    for Tf in (0.2, 0.5, 1, 2, 5, 10, 30, 100):
        # self-consistent mass: m = R * mp * sheet_cap / gap_cap(m)
        # solve m * cap(m) = R mp * sheet_cap; take the lowest root (m*cap has a maximum when Boltzmann suppression sets in)
        need = R_obs*mp*sheet_cap(Tf); mg = np.geomspace(0.05, 3000, 600)
        gv = np.array([mm*cap(mm, Tf, 2, qchi) for mm in mg]) - need
        idx = np.where(np.sign(gv[:-1]) != np.sign(gv[1:]))[0]
        m = brentq(lambda mm: mm*cap(mm, Tf, 2, qchi) - need, mg[idx[0]], mg[idx[0] + 1]) if len(idx) else np.nan
        H = 1.66*np.sqrt(gstar(Tf))*Tf**2/MPl
        Gc = np.sqrt(H/Tf**5)
        mu = m*mp/(m + mp)
        sig = Gc**2*mu**2/np.pi*0.3894e-27      # GeV^-2 -> cm^2
        rows.append((qchi, Tf, m, Gc, sig))
        out.append(f"   {Tf:8.1f}    {label:20s} {m:10.2f}           {Gc:10.2e}          {sig:10.1e}")
    out.append("")
mq1 = [r[2] for r in rows if r[0] == 1.0 and np.isfinite(r[2])]
# threshold: lowest T_f at which the simplest content can hold enough charge
def has_sol(Tf, q=1.0):
    need = R_obs*mp*sheet_cap(Tf); mg = np.geomspace(0.05, 3000, 400)
    return max(mm*cap(mm, Tf, 2, q) for mm in mg) >= need
Tmin = brentq(lambda T: (1 if has_sol(T) else -1) + 1e-9*T, 2, 10, xtol=0.05) if not has_sol(2) else 2
sig_hi = max(r[4] for r in rows if r[0] == 1.0 and np.isfinite(r[4]))
out += [f"SIMPLEST CONTENT (q = 1): Ocean particle mass {min(mq1):.1f}-{max(mq1):.1f} GeV (~{min(mq1)/mp:.0f}-{max(mq1)/mp:.0f} proton masses).",
        f"Crossing must shut ABOVE T_f ~ {Tmin:.1f} GeV; below that the Ocean particle is too heavy/slow to hold enough charge and",
        "the ratio 5.36 cannot be reached for ANY mass (q = 1/3 content needs T_f > ~100 GeV and m ~ 100 GeV).",
        "",
        "PRE-REGISTERED CLAIMS (written from the numbers above, BEFORE looking up detector limits):",
        f" P25a: the Ocean is made of particles of ~8.5-13 GeV (simplest content; 12.6 GeV at the threshold T_f = 7.3); q = 1/3 content: ~100 GeV.",
        " P25b: no annihilation signal from galaxy centres (no anti-Ocean survives).",
        f" P25c: contact coupling -> nucleon cross-section <= ~{sig_hi:.0e} cm^2 (order of magnitude): far below the neutrino",
        "       floor, so underground detectors should see NOTHING from the Ocean. (Weak test: only a detection would hurt.)",
        " Kill: a dark-matter detection at a clearly different mass (e.g. 1 TeV, or a keV line as all of DM) or above the",
        "       P25c strength at ~9 GeV with contact coupling; a convincing Galactic-centre annihilation signal.",
        " Not an explanation of 5.36: it is traded for the mass (~9 GeV). It would be explained only if the grid fixed the",
        "       Ocean particle's mass at ~9-10 proton masses."]
txt = "\n".join(out); print(txt); open("iter25_predictions.txt", "w").write(txt + "\n")
