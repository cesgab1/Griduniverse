"""
ITERATION 27: full chemical-equilibrium sharing for the twin Ocean (re-test of iteration 26 B, window FIXED beforehand at commit
7043860: twin b' mass 12.5-21 GeV, i.e. f/v = 3-5).
Conserved while crossing is open (T_f between 5 and 100 GeV, below the electroweak sphaleron era, so no SM B+L violation):
  B_tot = B_SM + B_twin   (crossing moves baryon number between sheets and gaps)
  L_SM                    (two choices, counted: L = 0, or the standard sphaleron leftover L = -(51/28) B_SM ... evaluated
                           self-consistently as L/B_tot fixed at -51/28)
  Q_SM = 0                (electric neutrality of our world)
  twin: no twin light force (iteration 26 A), so b' carries only B_twin.
Every species i: mu_i = b_i mu_B + l_i mu_L + q_i mu_Q; asymmetry dn_i = g_i c(m_i/T) mu_i T^2/6 (bosons: c -> 2 at m = 0).
Masses included (top, W, b, c, tau ... Boltzmann-suppressed automatically). m_DM = 3 m_b' (as iteration 26).
PASS if the solved m_b' lies in 12.5-21 GeV for some T_f in 5-100 GeV.
"""
import numpy as np, importlib.util, io, contextlib
from scipy.optimize import brentq
with contextlib.redirect_stdout(io.StringIO()):
    spec = importlib.util.spec_from_file_location("i25", "iter25_crossing.py"); i25 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i25)
cap, R, mp = i25.cap, i25.R_obs, i25.mp

# (name, mass, g, b, l, q, boson)
SM = [("u", .0022, 6, 1/3, 0, 2/3, 0), ("c", 1.27, 6, 1/3, 0, 2/3, 0), ("t", 172.7, 6, 1/3, 0, 2/3, 0),
      ("d", .0047, 6, 1/3, 0, -1/3, 0), ("s", .095, 6, 1/3, 0, -1/3, 0), ("b", 4.18, 6, 1/3, 0, -1/3, 0),
      ("e", .000511, 2, 0, 1, -1, 0), ("mu", .1057, 2, 0, 1, -1, 0), ("tau", 1.777, 2, 0, 1, -1, 0),
      ("nue", 0, 1, 0, 1, 0, 0), ("numu", 0, 1, 0, 1, 0, 0), ("nutau", 0, 1, 0, 1, 0, 0),
      ("W", 80.4, 3, 0, 0, 1, 1)]

def weights(T):
    """unit-mu asymmetry weights w_i = g_i c_i (cap() with q=1 gives g*c)"""
    return [cap(m, T, g, 1.0, boson=bool(bo)) for (_, m, g, b, l, q, bo) in SM]

def solve(T, mb, Lmode):
    w = weights(T)
    bvec = np.array([s[3] for s in SM]); lvec = np.array([s[4] for s in SM]); qvec = np.array([s[5] for s in SM]); w = np.array(w)
    # charges as linear functions of (mu_B, mu_L, mu_Q):  X = sum_i x_i w_i (b_i mu_B + l_i mu_L + q_i mu_Q)
    M = lambda x: np.array([np.sum(x*w*bvec), np.sum(x*w*lvec), np.sum(x*w*qvec)])
    Brow, Lrow, Qrow = M(bvec), M(lvec), M(qvec)
    wtw = cap(mb, T, 6, 1.0)                     # twin b': g = 6, b = 1/3
    Btw_row = np.array([wtw/9, 0, 0])            # B_twin = (1/3) * w * (mu_B/3)
    Btot_row = Brow + Btw_row
    # equations: Q = 0 ; L = 0  or  L = -(51/28) B_tot ; normalise mu_B = 1
    A = np.array([Qrow, Lrow if Lmode == "L0" else Lrow + 51/28*Btot_row])
    # unknowns mu_L, mu_Q with mu_B = 1
    sol = np.linalg.solve(A[:, 1:], -A[:, 0])
    mu = np.array([1.0, *sol])
    return (Btw_row @ mu)/(Brow @ mu)           # n_twin baryons / n_SM baryons

out = ["ITERATION 27: full sharing (neutrality, lepton number, masses, W) -- re-test of the twin weight window 12.5-21 GeV", "",
       "   T_f [GeV]   L choice     needed m_b' [GeV]   m_DM [GeV]   (crude iteration-26 value)"]
crude = {10: 9.54, 20: 8.62, 30: 8.56, 50: 8.84, 100: 9.50}
res = []
for Lmode in ("L0", "Lsph"):
    for T in (5, 7, 10, 20, 30, 50, 100):
        F = lambda mb: solve(T, mb, Lmode)*3*mb/mp - R
        mg = np.geomspace(0.1, 500, 300); gv = np.array([F(x) for x in mg])
        idx = np.where(np.sign(gv[:-1]) != np.sign(gv[1:]))[0]
        mb = brentq(F, mg[idx[0]], mg[idx[0] + 1]) if len(idx) else np.nan
        res.append((Lmode, T, mb))
        out.append(f"   {T:6.0f}      {Lmode:5s}       {mb:10.2f}        {3*mb:7.1f}       {crude.get(T, '')}")
    out.append("")
inwin = [(L, T, round(m, 2)) for L, T, m in res if np.isfinite(m) and 12.5 <= m <= 21]
out.append("RESULT: " + (f"inside the pre-set window for {inwin}" if inwin else "NEVER inside the pre-set window 12.5-21 GeV"))
fin = [m for _, _, m in res if np.isfinite(m)]
out.append(f"range of needed m_b' over all choices: {min(fin):.2f}-{max(fin):.2f} GeV  ->  f/v = {min(fin)/4.18:.2f}-{max(fin)/4.18:.2f}")
txt = "\n".join(out); print(txt); open("iter27_full_sharing.txt", "w").write(txt + "\n")
