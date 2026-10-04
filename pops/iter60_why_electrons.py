"""ITERATION 60 (pre-registered in PREREG_60.md): size vs own pop size for the charged particles present at payday."""
hbar_c = 197.3269804e-15   # MeV m
parts = [("electron", 0.51099895, 1e-18, "upper limit (point-like)"), ("proton", 938.272, 0.8409e-15, "charge radius"),
         ("helium-4 nucleus", 3727.379, 1.6755e-15, "charge radius"), ("deuteron", 1875.613, 2.12799e-15, "charge radius")]
out = ["ITERATION 60: who can hold a pop open? (size must be smaller than own pop size hbar/mc; expectations pre-registered)", ""]
for n, m, r, note in parts:
    pop = hbar_c/m
    out.append(f" {n:18s} pop size {pop:.2e} m   own size {r:.2e} m ({note})   size/pop = {r/pop:.1e}   -> {'HOLDS a pop (pays)' if r < pop else 'cannot (never pays)'}")
out += ["", "Result: only electrons pay; protons and nuclei are larger than their own pop size. The iteration-58 rule follows from R1.",
        "Payday is collective (the gas's thermal tie to the light), so bound and free electrons pay together at z ~ 124.",
        "Neutral particles (neutrinos, dark matter) hold no Coulomb pop and never pay."]
txt = "\n".join(out); print(txt); open("iter60_why_electrons.txt", "w").write(txt + "\n")
