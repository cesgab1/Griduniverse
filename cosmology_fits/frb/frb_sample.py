import json, glob, os, numpy as np
D="/home/claude/frb_repo/frb/data/"
rows=[]
for f in sorted(glob.glob(D+"FRBs/FRB*.json")):
    d=json.load(open(f)); name=d["FRB"]
    if "DM" not in d: continue
    hs=glob.glob(D+f"Galaxies/{name.replace('FRB','')}/*host*.json")
    if not hs: continue
    h=json.load(open(hs[0])); r=h.get("redshift",{})
    z=r.get("z_spec", r.get("z")); spec="z_spec" in r
    if z is None: continue
    dmism=d.get("DMISM",{}).get("value", np.nan)
    rows.append((name, float(z), spec, d["DM"]["value"], dmism))
print(len(rows), "FRBs with host redshift;", sum(r[2] for r in rows), "spectroscopic")
np.save("frb_rows.npy", np.array(rows, dtype=object), allow_pickle=True)
for r in rows[:5]: print(r)
zs=np.array([r[1] for r in rows]); print("z range", zs.min(), zs.max(), " missing DMISM:", sum(np.isnan(r[4]) for r in rows))
