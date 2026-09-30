"""Do galaxies whose rotation curves sag (fitted neighbour pull e) actually have stronger neighbour pulls
measured independently from maps of surrounding galaxies (Chae et al. 2021, Table 3, log e_N,env, 109 galaxies)?"""
import re, json, numpy as np
from scipy.stats import spearmanr, kendalltau
txt = open("/home/claude/chae.txt").read()
norm = lambda s: s.replace(" ", "").upper()
# Table 3: two columns per line
i0 = txt.index("Newtonian Environmental Field Strength"); t3 = txt[i0:i0+20000]
env = {}
for m in re.finditer(r"([A-Z][A-Za-z0-9\-]*(?: [A-Z0-9][A-Za-z0-9\-]*)?)\s+([−\-]\d\.\d{3}) ± (\d\.\d{3})\s+([−\-]\d\.\d{3}) ± (\d\.\d{3})", t3):
    name = norm(m.group(1)); f = lambda s: float(s.replace("−", "-"))
    env[name] = dict(maxc=f(m.group(2)), emax=float(m.group(3)), noc=f(m.group(4)))
print(f"Table 3 parsed: {len(env)} galaxies")
# Table 2: fitted e-tilde (value line: name flag x0 e- err ...)
t2 = {}
for m in re.finditer(r"^([A-Z][A-Za-z0-9\-]*(?: [A-Z0-9][A-Za-z0-9\-]*)?)\s+([PABC])\s+([−\-]\d+\.\d+)\s+([−\-]?\d\.\d{3})-", txt, re.M):
    t2[norm(m.group(1))] = dict(q=m.group(2), x03=float(m.group(3).replace("−", "-")), e=float(m.group(4).replace("−", "-")))
print(f"Table 2 parsed: {len(t2)} galaxies")
ours = {norm(x["name"]): x for x in json.load(open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/neighbour_test.json"))}
def key(n):  # SPARC file names: DDO064, NGC0100, UGCA281, ESO079-G014 ...
    return n
out = {}
for lab, src, fld in (("our fits", ours, "e"), ("Chae et al. fits (Table 2, quality P, x0,3 < -10.6)", t2, "e")):
    xs, ys, names = [], [], []
    for n, v in env.items():
        k = n.replace("NGC0", "NGC0")
        if k in src:
            s = src[k]
            if src is t2 and not (s["q"] == "P" and s["x03"] < -10.6): continue
            xs.append(v["maxc"]); ys.append(s[fld]); names.append(k)
    xs, ys = np.array(xs), np.array(ys)
    rho, p = spearmanr(xs, ys); tau, pt = kendalltau(xs, ys)
    med = np.median(xs); lo, hi = ys[xs <= med], ys[xs > med]
    print(f"\n{lab}: {len(xs)} galaxies matched")
    print(f"   rank correlation of fitted sag e with mapped neighbour pull: Spearman {rho:+.2f} (p = {p:.3f}), Kendall {tau:+.2f} (p = {pt:.3f})")
    print(f"   quieter half of environments: median e {np.median(lo):.3f};  busier half: median e {np.median(hi):.3f}")
    out[lab] = dict(n=len(xs), rho=rho, p=p, lo=float(np.median(lo)), hi=float(np.median(hi)))
# amplitude check: MOND external field from Newtonian e_N in deep regime e ~ sqrt(e_N)
eN = np.array([v["maxc"] for v in env.values()]); eN0 = np.array([v["noc"] for v in env.values()])
print(f"\nexpected size: mapped e_N median 10^{np.median(eN):.2f} (max clustering) -> MOND-scale e ~ sqrt = {np.sqrt(10**np.median(eN)):.3f};"
      f" no clustering 10^{np.median(eN0):.2f} -> {np.sqrt(10**np.median(eN0)):.3f}.  Our median fitted e: {np.median([x['e'] for x in ours.values()]):.3f}")
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/neighbour_env.json", "w"), indent=1)
