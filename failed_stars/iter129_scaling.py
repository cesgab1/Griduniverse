"""PREREG 129 (+addendum): do brown-dwarf disks follow the stars' dust-mass line?"""
import numpy as np, pandas as pd
rng = np.random.default_rng(129)
h=6.626e-34; k=1.381e-23; c=2.998e8; pc=3.086e16; Me=5.972e24
def mdust_earth(F_mJy, d_pc, T, nu=337e9, kappa_cgs=2.3*(337/230)**0.4):
    B = 2*h*nu**3/c**2/np.expm1(h*nu/(k*T))
    return F_mJy*1e-29*(d_pc*pc)**2/(kappa_cgs*0.1*B)/Me
out=[]; P=lambda s="": (print(s), out.append(s))

oph = pd.read_csv("data/testi2022_tableA1_L1688_VLMO.csv", comment="#")
oph["T_A"] = 20.0; oph["T_B"] = 25.0*(10**oph.logL)**0.25
lup = pd.read_csv("data/sanchis2020_lupus_BD_stars.csv", comment="#").dropna(subset=["Mstar_Msun"])
sets = {}
for conv in ("A","B"):
    d = oph.assign(Md=mdust_earth(oph.Ftot_mJy, oph.d_pc, oph["T_"+conv]))
    sets[f"L1688 conv {conv}"] = (d.Mstar_Msun.values, d.Md.values, d.Ftot_upper.values.astype(bool), 1.9, 0.4)
sets["Lupus conv A"] = (lup.Mstar_Msun.values, lup.Mdust_Mearth.values, lup.Mdust_upper.values.astype(bool), 1.73, 0.25)

def analyse(name, M, Md, up, slope, s_err):
    star = M > 0.1; bd = M <= 0.08
    x, y = np.log10(M), np.log10(Md)
    P(f"\n{name}: stars {star.sum()} (det {(star&~up).sum()}), brown dwarfs {bd.sum()} (det {(bd&~up).sum()})")
    for lab, b_fixed in (("published slope", slope), ("own slope", None)):
        sd = star & ~up
        if b_fixed is None:
            if sd.sum() < 3: P("   own slope: too few stars"); continue
            b, a = np.polyfit(x[sd], y[sd], 1)
        else:
            b = b_fixed; a = np.median(y[sd] - b*x[sd])
        res = y - (a + b*x)
        for treat, sel in (("detections only", bd & ~up), ("limits as values", bd)):
            r = res[sel]
            if len(r) == 0: continue
            boots = []
            for _ in range(4000):   # bootstrap BDs and star normalisation; slope drawn from its error if published
                bb = b if b_fixed is None else rng.normal(slope, s_err)
                si = rng.choice(np.where(sd)[0], sd.sum()); aa = np.median(y[si]-bb*x[si])
                bi = rng.choice(np.where(sel)[0], sel.sum())
                boots.append(np.median(y[bi]-(aa+bb*x[bi])))
            med, err = np.median(r), np.std(boots)
            brk = abs(med) > 0.3 and abs(med) > 2*err
            P(f"   {lab:15s} (b={b:.2f}) {treat:16s}: BD median residual {med:+.2f} +/- {err:.2f} dex "
              f"(n={len(r)}) -> {'BREAK-candidate' if brk else 'on the line'}")
for n, args in sets.items(): analyse(n, *args)
usco = pd.read_csv("data/vanderplas2016_USco_Oph_BD.csv", comment="#")
u = usco[usco.region=="UpperSco"]
P(f"\nUpper Sco (older, BDs only, descriptive): median dust {np.median(u.Mdust_A13):.2f} (A13 T) / "
  f"{np.median(u.Mdust_vdP16):.2f} (vdP16 T) Earth masses at median M = {np.median(u.Mstar_Msun):.2f} Msun")
open("RESULT_129_numbers.txt","w").write("\n".join(out)+"\n")
