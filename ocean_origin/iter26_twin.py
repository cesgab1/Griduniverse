"""
ITERATION 26 (pre-registration part): do the gaps hold a TWIN copy of matter (twin strong force, 'fraternal' twin with only the
third generation)?  BORROWED: Twin Higgs asymmetric DM (Garcia Garcia, Lasenby, March-Russell PRL 115 121801; Farina 1506.03520).
Grid reading: the sheets hold our matter; the gaps hold a twin copy whose pieces come out heavier. Kept kernel of Coalesce's
'anglerfish' idea: heavier because built from the same kind of pieces, bundled.

Checks and PASS/FAIL criteria (written before looking up the newest N_eff measurements):
 A. Dark radiation: twin neutrino (always) and twin photon (only if the twin has its own light force). After the two worlds stop
    trading heat at T_dec, the twin radiation is diluted by everything our world later annihilates.
    PASS if Delta N_eff is below the current 95% upper limit (looked up AFTER this commit). Also compared with our own
    pre-registered window for the jostling radiation (P3: < 0.107): the two share one budget.
 B. Consistency with iteration 25: the twin charge is carried by twin b' quarks (3 colours x 2 spins, charge 1/3) at the time
    sharing stops. Solve for the twin b' mass that gives omega_c/omega_b = 5.36, with m_DM = 3 m_b' (heavy-quark baryon;
    free choice). Twin Higgs expectation: m_b' = m_b * f/v with f/v = 3-5 -> 12.5-21 GeV.
    PASS if the solved m_b' lands in 12.5-21 GeV for some T_f in 1-100 GeV (loose test: T_f is free).
 C. Self-collisions: sigma ~ 4 pi / Lambda'^2 (strong-force size), m_DM from B, Lambda' = 1-5 GeV.
    PASS if sigma/m < 1 cm^2/g (Bullet Cluster family).
 D. Dark disks: only the version with a twin photon can cool into disks; it is flagged, not computed.
Free choices: fraternal twin content; m_DM = 3 m_b'; T_dec scanned 1-10 GeV; contact sharing as in iteration 25.
"""
import numpy as np, importlib.util, io, contextlib
from scipy.optimize import brentq
with contextlib.redirect_stdout(io.StringIO()):
    spec = importlib.util.spec_from_file_location("i25", "iter25_crossing.py"); i25 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i25)
cap, sheet_cap, R, mp = i25.cap, i25.sheet_cap, i25.R_obs, i25.mp

def gs_SM(T):   # entropy dof of our world (rough step function, standard values)
    for Tc, g in ((0.0005, 3.909), (0.15, 10.75), (1.0, 61.75), (4.0, 75.75), (80, 86.25)):
        if T < Tc: return g
    return 106.75
out = ["ITERATION 26: twin copy of matter in the gaps -- checks (criteria committed before the newest N_eff lookup)", "",
       "A. Dark radiation Delta N_eff (twin tau' heavy at T_dec; twin glueballs decay into OUR world, not twin radiation)",
       "   T_dec [GeV]   g_*s(ours)   T'/T today   neutrino-only   with twin photon"]
dN = {}
for Td in (1.0, 2.0, 3.0, 5.0, 10.0):
    r = (3.909/gs_SM(Td))**(1/3)                  # twin light sector entropy conserved separately; ours dumps into photons
    f = (11/4)**(4/3)*r**4*4/7
    dN[Td] = (f*1.75, f*3.75)
    out.append(f"   {Td:8.1f}      {gs_SM(Td):6.2f}       {r:.3f}        {f*1.75:.3f}          {f*3.75:.3f}")
out += ["   (T_dec below the QCD transition, ~0.15 GeV, would give 0.5-1: excluded outright; twin worlds need early decoupling)", ""]

out += ["B. Twin b' mass needed for omega_c/omega_b = 5.36 (m_DM = 3 m_b'); twin Higgs expects 12.5-21 GeV",
        "   T_f [GeV]   needed m_b' [GeV]   m_DM [GeV]"]
sols = []
for Tf in (1, 2, 3, 5, 7, 10, 20, 30, 50, 100):
    need = R*mp*sheet_cap(Tf)
    F = lambda mb: (3*mb)*cap(mb, Tf, 6, 1/3) - need
    mg = np.geomspace(0.05, 2000, 500); gv = np.array([F(x) for x in mg])
    idx = np.where(np.sign(gv[:-1]) != np.sign(gv[1:]))[0]
    mb = brentq(F, mg[idx[0]], mg[idx[0] + 1]) if len(idx) else np.nan
    sols.append((Tf, mb)); out.append(f"   {Tf:6.0f}      {mb:10.2f}          {3*mb:8.1f}")
ok = [(t, m) for t, m in sols if np.isfinite(m) and 12.5 <= m <= 21]
out += [f"   -> in the twin-Higgs window for T_f = {[t for t, _ in ok]} GeV" if ok else "   -> never in the twin-Higgs window", ""]

out += ["C. Self-collisions sigma/m (sigma ~ 4 pi / Lambda'^2)"]
GeV2cm2, GeVg = 0.3894e-27, 1.783e-24
for mdm in (30.0, 60.0):
    for Lp in (1.0, 2.0, 5.0):
        s = 4*np.pi/Lp**2*GeV2cm2
        out.append(f"   m_DM {mdm:4.0f} GeV, Lambda' {Lp:.0f} GeV: sigma/m = {s/(mdm*GeVg):.1e} cm^2/g")
out += ["", "D. Twin photon version: twin baryons are charged -> twin atoms can cool -> dark disks (Gaia limits ~ few % of DM in a",
        "   thin disk). Flagged; the neutrino-only (no twin light force) version avoids it.", ""]
txt = "\n".join(out); print(txt); open("iter26_part1.txt", "w").write(txt + "\n")
