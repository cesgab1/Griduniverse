"""Joint profile of counting exponent beta and link-stretch exponent s; law d ln rho_DE/d ln a = beta (q + 1 - s)."""
import json, numpy as np, collections
from scipy.interpolate import RectBivariateSpline
R = collections.defaultdict(dict); L = {}
import os
for line in open("joint_raw.txt").readlines() + (open("joint_raw_ext.txt").readlines() if os.path.exists("joint_raw_ext.txt") else []):
    parts = line.split(" ", 3)
    if len(parts) < 4 or "ONLYRESULT" not in parts[3]: continue
    sn, b, s, rest = parts; d = json.loads(rest.split("ONLYRESULT ", 1)[1])
    if b == "LCDM": L[sn] = d["LCDM"]["chi2"]
    else: R[sn][(float(b), float(s))] = d["LAW"]["chi2"]
out = ["Joint fit of beta (counting exponent; sqrt(N) -> 1/2) and s (link-stretch exponent); simplified q; Delta chi2 vs LCDM", ""]
for sn in ("PANTHEON", "DESY5", "UNION3"):
    B = sorted({k[0] for k in R[sn]}); S = sorted({k[1] for k in R[sn]})
    Z = np.array([[R[sn].get((b, s), np.nan) - L[sn] for s in S] for b in B])
    out.append(f"{sn}: rows beta {B}, columns s {S}")
    for b, row in zip(B, Z): out.append(f"  beta {b:5.3f}: " + " ".join(f"{x:+6.1f}" for x in row))
    sp = RectBivariateSpline(B, S, Z, kx=3, ky=3)
    bb = np.linspace(B[0], B[-1], 301); ss = np.linspace(S[0], S[-1], 301); G = sp(bb, ss)
    i, j = np.unravel_index(np.argmin(G), G.shape); gmin = G[i, j]
    inside = G - gmin < 2.30                                   # 68% region, 2 parameters
    out.append(f"  best: beta = {bb[i]:.2f}, s = {ss[j]:.2f}, Delta chi2 vs LCDM {gmin:+.1f};  68% region: beta {bb[inside.any(1)].min():.2f}-{bb[inside.any(1)].max():.2f}, s {ss[inside.any(0)].min():.2f}-{ss[inside.any(0)].max():.2f}")
    k = np.argmin(np.abs(bb - 0.5)); l1 = np.argmin(np.abs(ss - 1.0))
    out.append(f"  Claim 1 point (beta 1/2, s 1): {G[k, l1] - gmin:+.2f} above the best (2-parameter 68% = 2.30, 95% = 6.18)")
    # degeneracy direction: the quantity the data fix
    pts = np.argwhere(inside); bv, sv = bb[pts[:, 0]], ss[pts[:, 1]]
    c = np.corrcoef(bv, sv)[0, 1]; out.append(f"  beta-s correlation inside 68% region: {c:+.2f}")
    # what data fix: today's slope beta(q0+1-s) and matter-era slope beta(1.5-s)
    out.append(f"  inside 68%: today's slope beta(1 + q0 - s) with q0 = -0.53: {np.min(bv*(0.47 - sv)):.2f} to {np.max(bv*(0.47 - sv)):.2f};"
               f"  matter-era slope beta(1.5 - s): {np.min(bv*(1.5 - sv)):.2f} to {np.max(bv*(1.5 - sv)):.2f}  (Claim 1: {0.5*(0.47-1):.2f}, {0.5*0.5:.2f})")
    out.append("")
txt = "\n".join(out); print(txt); open("joint_scan.txt", "w").write(txt + "\n")
