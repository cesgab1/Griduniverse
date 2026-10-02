"""
Figure 10: one cell of the grid, and how big cells may be.
(a) A real cell of a random 3-D mosaic (Poisson-Voronoi), with its links: one link through each shared face to a neighbour.
(b) Cells are not all the same size: volume and face-count spread over ~3000 cells (computed).
(c) The allowed size window: cell size must be below 5.7e-28 m (light shows no grid effect: LHAASO GRB 221009A, our
    lorentz/lorentz_check.txt) and, in the model, near the Planck length (the density cap that turns collapse into a
    bounce sits at ~0.4 Planck density).
Outputs fig10_grid_cell.png/.pdf (dark style to match the cover).
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, ConvexHull
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(3)
BG, GRID, GOLD, LIGHT, VIO, TXT, MUT = "black", "#58aaff", "#ffc94d", "#c4f1ff", "#b48cff", "#e6ecf7", "#93a2bd"
plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": TXT, "axes.labelcolor": TXT, "xtick.color": MUT,
                     "ytick.color": MUT, "axes.edgecolor": "#3a4660", "font.size": 10})
AZ, EL = np.radians(35), np.radians(22)
def proj(p):
    x, y, z = p[..., 0], p[..., 1], p[..., 2]
    xr = x*np.cos(AZ) - y*np.sin(AZ); yr = x*np.sin(AZ) + y*np.cos(AZ)
    d = yr*np.cos(EL) - z*np.sin(EL); up = z*np.cos(EL) + yr*np.sin(EL); s = 6/(6 + d)
    return np.stack([xr*s, up*s], -1), d
def glow(ax, P2, color, lw=1.0, alpha=1.0, zo=3, ls="-"):
    for k, a in ((8, .035), (4, .09), (2, .25), (1, 1)):
        ax.plot(P2[:, 0], P2[:, 1], ls, color=color, lw=lw*k, alpha=a*alpha, zorder=zo, solid_capstyle="round")

# random 3-D mosaic: 4000 cells in a periodic box of side 16 (mean cell volume ~1)
L = 16.0; n = 4096; P = rng.uniform(0, L, (n, 3))
shifts = np.array([(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])*L
Q = np.vstack([P + s for s in shifts]); home = 13*n                     # index offset of the unshifted copy
vor = Voronoi(Q)
vols, faces = [], []
for i in range(n):
    reg = vor.regions[vor.point_region[home + i]]
    if -1 in reg or len(reg) == 0: continue
    V = vor.vertices[reg]; vols.append(ConvexHull(V).volume)
vols = np.array(vols); mv = vols.mean(); vols /= mv
# face counts: ridges per cell
cnt = np.zeros(len(Q), int)
for a, b in vor.ridge_points: cnt[a] += 1; cnt[b] += 1
faces = cnt[home:home + n]

fig = plt.figure(figsize=(14, 5.6), facecolor=BG)
# (a) one cell
ax = fig.add_axes([0.0, 0.1, 0.36, 0.78]); ax.set_facecolor(BG); ax.axis("off")
c0 = home + int(np.argmin(np.linalg.norm(P - L/2, axis=1))); centre = Q[c0]
ridges = [(rv, rp) for rv, rp in zip(vor.ridge_vertices, vor.ridge_points) if c0 in rp and -1 not in rv]
allv = np.vstack([vor.vertices[rv] for rv, _ in ridges]); scale = 1/np.max(np.linalg.norm(allv - centre, axis=1))
def T(p): return (p - centre)*scale
# faces (back first): order polygon vertices, draw edges with glow; translucent fill
from matplotlib.patches import Polygon
polys = []
for rv, rp in ridges:
    V = T(vor.vertices[rv]); c = V.mean(0); nrm = np.cross(V[1] - V[0], V[2] - V[0]); nrm /= np.linalg.norm(nrm) + 1e-12
    u = V[0] - c; u /= np.linalg.norm(u); w = np.cross(nrm, u); ang = np.arctan2((V - c) @ w, (V - c) @ u)
    V = V[np.argsort(ang)]; P2, d = proj(np.vstack([V, V[:1]])); polys.append((d.mean(), P2, rp))
polys.sort(key=lambda t: -t[0])
for dm, P2, rp in polys:
    ax.add_patch(Polygon(P2[:-1], closed=True, fc=GRID, alpha=0.05, ec="none", zorder=2))
    glow(ax, P2, GRID, lw=0.9, alpha=0.55 if dm > 0 else 1.0, zo=3)
# links to neighbours (centre-to-centre through each face)
for dm, P2, rp in polys:
    nb = rp[1] if rp[0] == c0 else rp[0]; tn = T(Q[nb]); seg = np.vstack([np.zeros(3), tn*min(1.0, 1.25/np.linalg.norm(tn))]); S2, _ = proj(seg)
    glow(ax, S2, GOLD, lw=0.8, alpha=0.5 if dm > 0 else 0.95, zo=4)
    ax.scatter(*S2[1], s=16, color=GOLD, zorder=5, lw=0)
o2, _ = proj(np.zeros((1, 3)))
for s, a in ((600, .06), (220, .2), (60, 1)): ax.scatter(*o2.T, s=s, color="#ffffff" if s == 60 else VIO, alpha=a, lw=0, zorder=6)
ax.set_xlim(-1.4, 1.4); ax.set_ylim(-1.3, 1.25); ax.set_aspect("equal")
ax.text(-1.37, 1.22, "(a) One cell of the grid (a real random cell)", fontsize=12, fontweight="bold", va="top")
ax.text(-1.37, 1.08, f"{len(ridges)} faces; average over all cells {faces.mean():.1f} (a cube has 6)", fontsize=9, color=MUT, va="top")
ax.text(-1.37, -1.27, "GRID CELL (Planck cell). White: its NODE. Blue: its WALLS, each shared with one\nneighbouring cell. Gold: its TENSION LINKS, one through each wall to the\nneighbouring node; a link's strength is set by its wall's area.", fontsize=9, color=TXT, va="top")
# (b) size spread
ax = fig.add_axes([0.42, 0.14, 0.25, 0.66]); ax.set_facecolor(BG)
ax.hist(vols, bins=40, range=(0, 2.6), color=GRID, alpha=0.85, edgecolor=BG)
ax.set_xlabel("cell volume / average"); ax.set_ylabel("number of cells")
for s in ("top", "right"): ax.spines[s].set_visible(False)
p1, p99 = np.percentile(vols, [1, 99])
ax.set_title("(b) Cells come in a spread of sizes", loc="left", fontsize=12, fontweight="bold", color=TXT, pad=14)
ax.text(0.98, 0.95, f"{len(vols)} cells (computed)\nspread ±{vols.std():.0%} of average\n98% between {p1:.2f}x and {p99:.2f}x\nfaces per cell: {faces.min()}-{faces.max()}",
        transform=ax.transAxes, ha="right", va="top", fontsize=9, color=TXT)
ax.text(0, -0.3, "No cell needs a set size or shape: only the AVERAGE size\nmatters, and randomness is what keeps every direction equal.", transform=ax.transAxes, fontsize=9, color=MUT, va="top")
# (c) size window
ax = fig.add_axes([0.71, 0.14, 0.27, 0.66]); ax.set_facecolor(BG)
ax.set_xscale("log"); ax.set_xlim(1e-38, 1e-17); ax.set_ylim(0, 1); ax.set_yticks([])
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
lP, lmax = 1.616e-35, 5.7e-28
ax.axvspan(lP, lmax, color="#1d4f2f", alpha=0.9); ax.axvspan(lmax, 1e-17, color="#4a1c1c", alpha=0.9); ax.axvspan(1e-38, lP, color="#26262e", alpha=0.9)
ax.text(np.sqrt(lP*lmax), 0.88, "allowed\nwindow\n(7.5 powers\nof ten)", ha="center", va="top", fontsize=9, color="#bff0cf")
ax.text(3e-23, 0.9, "EXCLUDED:\nlight would show\nthe grid (colours\narrive at different\ntimes from far\ngamma-ray bursts)", ha="center", va="top", fontsize=8.5, color="#ffc9c9")
ax.text(1.6e-37, 0.9, "too\nsmall\nto\nmean\nany-\nthing", ha="center", va="top", fontsize=8.5, color=MUT)
for x, lab, y in ((lP, "MIN: Planck\nlength 1.6e-35 m\n(the model's\ncell size, set\nby its density\ncap)", 0.30), (lmax, "MAX allowed\n5.7e-28 m", 0.30),
                  (1e-19, "LHC sees\ndown to\n~1e-19 m", 0.55)):
    ax.axvline(x, color=TXT if x in (lP, lmax) else MUT, lw=1.3 if x in (lP, lmax) else 0.8, ls="-" if x in (lP, lmax) else ":")
    ax.text(x*1.6, y, lab, fontsize=8, color=TXT if x in (lP, lmax) else MUT, va="top")
ax.set_xlabel("size of one cell (metres)")
ax.set_title("(c) How big can a cell be?", loc="left", fontsize=12, fontweight="bold", color=TXT, pad=14)
fig.text(0.005, 0.975, "Figure 10. The grid's building block: a random cell, its links, and the size it must have",
         fontsize=13, fontweight="bold", color=TXT, va="top")
for ext in ("png", "pdf"): fig.savefig(os.path.join(HERE, f"fig10_grid_cell.{ext}"), dpi=170, facecolor=BG)
print(f"cells {len(vols)}, spread {vols.std():.3f}, faces mean {faces.mean():.2f}, this cell {len(ridges)}")
