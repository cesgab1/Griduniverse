import numpy as np, os
from astropy.coordinates import SkyCoord; import astropy.units as u
here = os.path.dirname(os.path.abspath(__file__))
E = [("CMB",264.0,48.3,"v"),("RADIO",248,44,"v"),("QSO",238.2,28.8,"v"),("LOWL",239,64.3,"a"),("HPA",227,-15,"v"),
     ("PARITY",264,-17,"a"),("CLUSTER",280,-15,"v"),("BULK",282,6,"v"),("SN",316.5,12.6,"v"),("PARAM",48.8,-5.6,"v"),("COLD",209,-57,"v")]
def uv(l, b):
    l, b = np.radians(l), np.radians(b); return np.array([np.cos(b)*np.cos(l), np.cos(b)*np.sin(l), np.sin(b)])
V = np.array([uv(l, b) for _, l, b, _ in E]); axis = np.array([t == "a" for *_, t in E])
def stat(V, axis_mask):
    # iterate: choose signs of free entries to best align with current mean
    s = V.copy(); m = s.sum(0)
    for _ in range(20):
        dots = s @ m; flip = axis_mask & (dots < 0); s[flip] *= -1; m = s.sum(0)
    # try starting from each entry to avoid local optimum
    best = np.linalg.norm(m)/len(V); bestm = m
    for k in range(len(V)):
        s = V.copy(); m = s[k].copy()
        for _ in range(20):
            dots = s @ m; s[axis_mask & (dots < 0)] *= -1; m = s.sum(0)
        if np.linalg.norm(m)/len(V) > best: best, bestm = np.linalg.norm(m)/len(V), m
    return best, bestm/np.linalg.norm(bestm)
rng = np.random.default_rng(7); out = []
for ver, mask in (("Version 1 (vectors as published, axes free)", axis), ("Version 2 (all sign-free)", np.ones(len(E), bool))):
    R, d = stat(V, mask)
    sims = []
    for _ in range(20000):
        X = rng.normal(size=(len(E), 3)); X /= np.linalg.norm(X, axis=1)[:, None]; sims.append(stat(X, mask)[0])
    p = np.mean(np.array(sims) >= R)
    lc, bc = np.degrees(np.arctan2(d[1], d[0])) % 360, np.degrees(np.arcsin(d[2]))
    ang = [np.degrees(np.arccos(min(1, abs(v @ d) if m_ else v @ d))) for v, m_ in zip(V, mask)]
    out.append(f"{ver}: R = {R:.3f}, chance p = {p:.4f}, common direction (l,b) = ({lc:.0f}, {bc:+.0f})")
    out.append("   angle of each to the common direction: " + ", ".join(f"{E[i][0]} {a:.0f}" for i, a in enumerate(ang)))
    out.append(f"   within 45 deg: {sum(a < 45 for a in ang)} of {len(E)}")
cmb = uv(264.0, 48.3)
out.append("angles to the CMB dipole: " + ", ".join(f"{n} {np.degrees(np.arccos(np.clip(uv(l,b)@cmb,-1,1))):.0f}" for n, l, b, _ in E))
txt = "\n".join(out); print(txt); open(os.path.join(here, "iter104_align.txt"), "w").write(txt + "\n")
