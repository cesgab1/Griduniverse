import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
d = json.load(open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/crack_grid.json"))
fig, axs = plt.subplots(1, 3, figsize=(13.5, 5.4), facecolor="#fcfcfb")
titles = {"(a) horizontal/vertical cracks": "Horizontal/vertical cracks (your picture)", "(b) random-angle cracks": "Same cracking, random angles",
          "(c) random-point grid": "Random-point grid (framework's grid)"}
for ax, (lab, g) in zip(axs, d.items()):
    P = np.array(g["P"]); E = np.array(g["edges"]); src = g["src"]; D = {int(k): v for k, v in g["D"].items()}
    ax.add_collection(LineCollection(P[E], colors="#b9b8b2", linewidths=0.45))
    ids = np.array(list(D)); t = np.array([D[i] for i in ids]); r = np.linalg.norm(P[ids]-P[src], axis=1)
    T = np.median(t[(r > 0.27) & (r < 0.33)]); pts = []
    for u, v in E:
        if u in D and v in D and (D[u]-T)*(D[v]-T) < 0:
            f = (T-D[u])/(D[v]-D[u]); pts.append(P[u] + f*(P[v]-P[u]))
    pts = np.array(pts); circ = np.linspace(0, 2*np.pi, 300); rr = np.median(np.linalg.norm(pts-P[src], axis=1))
    ax.scatter(pts[:, 0], pts[:, 1], s=9, color="#2a78d6", zorder=3, linewidths=0)
    ax.plot(P[src][0]+rr*np.cos(circ), P[src][1]+rr*np.sin(circ), color="#52514e", lw=1, ls=(0, (3, 3)))
    ax.plot(*P[src], "o", color="#0b0b0b", ms=5, zorder=4)
    ax.set_xlim(0.05, 0.95); ax.set_ylim(0.05, 0.95); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for s_ in ax.spines.values(): s_.set_visible(False)
    ax.set_title(titles[lab], fontsize=11.5, color="#0b0b0b", pad=10)
fig.text(0.5, 0.03, "Blue dots: where a signal sent from the black dot along the links has reached at one moment. Dashed: a perfect circle for comparison.",
         ha="center", fontsize=10, color="#52514e")
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.08, wspace=0.04)
plt.savefig("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/crack_grain.png", dpi=150, facecolor="#fcfcfb")
