"""ITERATION 59 (pre-registered in PREREG_59.md): grid quantities vs the Ocean particle mass window 8.5-13 GeV."""
import numpy as np
mp, me, a = 0.938272, 0.000510999, 1/137.035999
cands = [(f"m_p x 2^{k} (layers)", mp*2**k) for k in range(1, 6)]
cands += [("m_p / sqrt(alpha)", mp/np.sqrt(a)), ("m_p / alpha", mp/a)]
cands += [(f"m_e / alpha^{n}", me/a**n) for n in (2, 3, 4)]
cands += [("5.36 m_p (equal numbers)", 5.364*mp)]
out = ["ITERATION 59: Ocean particle mass from grid quantities (window 8.5-13 GeV; expectations pre-registered)", ""]
for n, v in cands:
    out.append(f" {n:28s} {v:9.3f} GeV   {'HIT' if 8.5 <= v <= 13 else ''}")
out += ["", "GR level: only 'cold and collisionless' -- passed (iteration 25 checklist). No number to compare."]
txt = "\n".join(out); print(txt); open("iter59_dm_two_levels.txt", "w").write(txt + "\n")
