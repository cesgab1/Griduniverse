"""Iteration 75: 'starts big, fades fast, then freezes' vs the early record. See PREREG_75.md."""
import numpy as np
from scipy.optimize import brentq
exec(open("iter74_bang_dial.py").read().split("S_today = ")[0])   # reuse model + g* steps
def T_of_a(a):
    T = T0K*kB_GeV/a
    for _ in range(5): T = T0K*kB_GeV/a*(3.91/gstar(T)[1])**(1/3)
    return T
aP = a_of_T(1.2209e19)
EP, shP, _ = E_and_share(aP); rtotP = EP**2
events = {"Planck": aP, "electroweak": a_of_T(160), "QCD": a_of_T(0.15), "e+e- annihilation": a_of_T(5e-4),
          "first 3 minutes": a_of_T(7e-5), "equality": brentq(lambda a: Om/a**3 - Orad(a), 1e-6, 1e-2), "transparency": 1/1091}
aBBN, aRec = events["first 3 minutes"], events["transparency"]
out = ["Iteration 75: starts at share 1/2 at the Planck time, fades as a^-n, freezes into the late-time value. PREREG_75.md", ""]
for n in (3, 4, 6):
    for endpoint in ("A constant", "B fading law"):
        start = 0.5*rtotP
        fade = lambda a: start*(a/aP)**(-n)
        late = (lambda a: ODE) if endpoint.startswith("A") else (lambda a: E_and_share(a)[2])
        g = lambda la: np.log(fade(np.exp(la))) - np.log(late(np.exp(la)))
        lo, hi = np.log(aP)+1e-9, 0.0
        if g(lo)*g(hi) > 0:
            out.append(f"n = {n}, {endpoint}: never reaches the late-time value before today -> IMPOSSIBLE"); continue
        a_s = np.exp(brentq(g, lo, hi))
        def share(a):
            E2 = E_and_share(a)[0]**2 - E_and_share(a)[2]
            r = fade(a) if a < a_s else late(a)
            return r/(E2 + r)
        sB, sR = share(aBBN), share(aRec)
        ok = sB <= 0.04 and sR <= 0.0036
        Ts = T_of_a(a_s)
        near = [k for k, v in events.items() if abs(np.log10(T_of_a(v)/Ts)) < np.log10(3)]
        out.append(f"n = {n}, {endpoint}: freezes at a = {a_s:.2e} (T = {Ts:.2e} GeV) | share at 3 min {sB:.1e} (limit 0.04), "
                   f"at transparency {sR:.1e} (limit 0.0036) -> {'PASSES' if ok else 'EXCLUDED'} | near event: {near or 'none'}")
out.append("")
out.append("Event temperatures (GeV): " + ", ".join(f"{k} {T_of_a(v):.2e}" for k, v in events.items()))
txt = "\n".join(out); print(txt); open("iter75_fade_freeze.txt", "w").write(txt + "\n")
