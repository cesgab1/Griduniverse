"""Pre-registration figure: predicted w(z) for Claim 1 and the toy (memory kappa = 3), with LCDM, and binned values."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = json.load(open("predictions.json"))["DESY5"]
C1, TOYC, GRAY, INK, MUT = "#2a78d6", "#eb6834", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#c9c8c2", "text.color": INK,
                     "axes.labelcolor": MUT, "xtick.color": MUT, "ytick.color": MUT})
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2), facecolor="#fcfcfb")
ax = axs[0]; ax.set_facecolor("#fcfcfb")
for name, col, lab in (("CLAIM1", C1, "Claim 1 (instant law)"), ("TOY", TOYC, "Memory version (κ = 3)")):
    z = np.array(R[name]["z"]); w = np.array(R[name]["w"]); ax.plot(z, w, color=col, lw=2)
    c = R[name]["crossing"][0]; ax.plot([c], [-1], "o", ms=8, color=col, mec="#fcfcfb", mew=2, zorder=5)
    ax.annotate(f"{lab}\ncrosses w = −1 at z = {c:.2f}", (c, -1), xytext=(c + 0.25, -0.93 if name == "CLAIM1" else -1.09),
                color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=MUT, lw=0.8))
ax.axhline(-1, color=GRAY, lw=1.5, ls="--"); ax.text(2.55, -0.995, "Λ (w = −1, never crosses)", color=MUT, fontsize=9, va="bottom", ha="right")
ax.set_xlim(0, 2.6); ax.set_ylim(-1.17, -0.86); ax.set_xlabel("redshift z"); ax.set_ylabel("dark-energy equation of state w")
ax.grid(axis="y", color="#ecebe6", lw=0.8); [ax.spines[s].set_visible(False) for s in ("top", "right")]
ax.set_title("Predicted w(z)", loc="left", fontsize=12, fontweight="bold", color=INK)
ax = axs[1]; ax.set_facecolor("#fcfcfb"); ax.set_axisbelow(True)
bins = list(R["CLAIM1"]["w_bins"].keys()); x = np.arange(len(bins)); wdt = 0.36
for i, (name, col) in enumerate((("CLAIM1", C1), ("TOY", TOYC))):
    v = np.array([R[name]["w_bins"][b] for b in bins]) + 1
    ax.bar(x + (i - 0.5)*wdt*1.06, v, wdt, color=col, label=("Claim 1" if name == "CLAIM1" else "Memory version (κ = 3)"))
ax.axhline(0, color=GRAY, lw=1.2); ax.set_xticks(x); ax.set_xticklabels(bins); ax.set_xlabel("redshift bin")
ax.set_ylabel("w + 1   (0 = Λ)"); ax.grid(axis="y", color="#ecebe6", lw=0.8); [ax.spines[s].set_visible(False) for s in ("top", "right")]
ax.legend(frameon=False, fontsize=9, loc="lower left"); ax.set_title("Predicted bin averages (what DESI DR3 / Euclid measure)", loc="left", fontsize=12, fontweight="bold", color=INK)
fig.text(0.01, -0.03, "Frozen 2026-10-03 before DESI DR3. Best fits to DESI DR2 BAO + Planck 2018 distance priors + DES-Dovekie (Pantheon+ and Union3 agree to ±0.002).",
         fontsize=8.5, color=MUT)
fig.tight_layout(); fig.savefig("prediction_wz.png", dpi=170, facecolor="#fcfcfb", bbox_inches="tight"); print("ok")
