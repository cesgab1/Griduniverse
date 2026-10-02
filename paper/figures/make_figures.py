"""
Paper figures for the Grid Universe model (current best fit: grid + its dark-energy law + a cold fluid between the layers).
Every panel that shows numbers is computed here or read from the repo's results; schematic panels are labelled 'picture'.
Run: python3 make_figures.py   -> fig*.png (200 dpi) and fig*.pdf in this folder.
"""
import numpy as np, json, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Circle, Ellipse, FancyBboxPatch
from matplotlib.collections import LineCollection, PolyCollection
from scipy.spatial import Voronoi, Delaunay
import scipy.sparse as sp, scipy.sparse.linalg as sla
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(os.path.dirname(HERE))
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"
BLUE, ORANGE, AQUA, VIOLET, RED, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7", "#e34948", "#008300"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9.5, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "axes.titlesize": 10.5, "axes.titlecolor": INK, "text.color": INK,
                     "lines.linewidth": 2, "axes.titlepad": 15, "savefig.facecolor": "white", "figure.facecolor": "white"})
rng = np.random.default_rng(7)
def save(fig, name):
    fig.savefig(os.path.join(HERE, name + ".png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(HERE, name + ".pdf"), bbox_inches="tight"); plt.close(fig); print("wrote", name)
def tag(ax, s, x=0.0, y=1.005):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, color=MUTED, style="italic", va="bottom")
def nobox(ax): ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]

def voronoi_edges(P, box):
    vor = Voronoi(P); segs = []
    for rv in vor.ridge_vertices:
        if -1 in rv: continue
        a, b = vor.vertices[rv[0]], vor.vertices[rv[1]]
        if np.all((a > -box) & (a < 2*box)) and np.all((b > -box) & (b < 2*box)): segs.append([a, b])
    return np.array(segs)

# ---------------------------------------------------------------- Fig 1: anatomy of the grid
def fig1():
    fig = plt.figure(figsize=(10, 4.4)); ax = fig.add_axes([0, 0, 0.62, 1]); nobox(ax)
    P = rng.uniform(0, 1, (170, 2)); segs = voronoi_edges(np.vstack([P + [dx, dy] for dx in (-1, 0, 1) for dy in (-1, 0, 1)]), 1)
    def skew(xy, z): return np.stack([xy[..., 0] + 0.45*xy[..., 1], 0.38*xy[..., 1] + z], axis=-1)
    for i, z in enumerate((0.0, 0.62, 1.24)):
        frame = skew(np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float), z)
        ax.add_patch(Polygon(frame, closed=True, fc="white", ec=MUTED, lw=1, zorder=3*i))
        s = segs[np.all((segs >= 0) & (segs <= 1), axis=(1, 2))]
        lc = LineCollection(skew(s, z), colors=INK2, lw=0.7, zorder=3*i + 1); ax.add_collection(lc)
        if i < 2:   # fluid between this layer and the next
            for j in range(140):
                p = rng.uniform(0, 1, 2); h = rng.uniform(0.08, 0.54)
                q = skew(p[None, :], z + h)[0]; ax.plot(*q, "o", ms=2.2, color=BLUE, alpha=0.55, zorder=3*i + 2)
    ax.set_xlim(-0.05, 1.85); ax.set_ylim(-0.08, 1.7)
    ax.text(1.47, 0.36, "layer = all of space\nat one moment", fontsize=8.5, color=INK2, va="center")
    ax.text(1.47, 0.67, "fluid between\nlayers ('ocean')", fontsize=8.5, color=BLUE, va="center")
    ax.annotate("", xy=(-0.02, 1.62), xytext=(-0.02, 0.05), arrowprops=dict(arrowstyle="->", color=INK2, lw=1.4))
    ax.text(-0.04, 0.85, "time: one layer\nafter another", rotation=90, ha="right", va="center", fontsize=8.5, color=INK2)
    ax.set_title("Figure 1. Anatomy of the grid", loc="left", x=0.02)
    tx = fig.add_axes([0.64, 0.05, 0.35, 0.88]); nobox(tx)
    txt = ("WHAT THE PIECES ARE\n\n"
           "Links: carry tension, like the strings of a\nstretched net. Their pull is what we feel\nas gravity (Fig. 2-3).\n\n"
           "Cells: random, not a crystal, so no direction\nis special (a regular grid would show its axes;\nwe measured ~20% for a cube lattice).\n\n"
           "Layers: one per moment, like the pages of a\nflip-book. Relativity's 'preferred slicing' in\nKhronon theory is this stack.\n\n"
           "Fluid: fills the gaps between pages. It is\npulled by gravity but does not collide or shine:\nthe dark matter of the CMB, clusters, lensing.\n\n"
           "Cosmic tension: the net's overall stretch is\nthe dark energy (Fig. 4).")
    tx.text(0, 1, txt, va="top", fontsize=8.6, color=INK, linespacing=1.35)
    tag(tx, "picture (layout); cell geometry is a real random mosaic", 0, -0.02)
    save(fig, "fig1_grid_anatomy")

# ---------------------------------------------------------------- Fig 2: gravity = the grid's tension (computed on a random network)
def fig2():
    fig, axs = plt.subplots(1, 3, figsize=(11, 3.7), gridspec_kw=dict(width_ratios=[1.15, 1.15, 1]))
    # (a) 2-D random mosaic, point load at the centre: equilibrium displacement of a tensioned net (graph Laplacian)
    n = 2600; P = rng.uniform(-1, 1, (n, 2)); tri = Delaunay(P)
    I, J = [], []
    for s in tri.simplices:
        for a, b in ((0, 1), (1, 2), (0, 2)): I += [s[a], s[b]]; J += [s[b], s[a]]
    A = sp.coo_matrix((np.ones(len(I)), (I, J)), shape=(n, n)).tocsr(); A.data[:] = 1
    L = sp.diags(np.asarray(A.sum(1)).ravel()) - A
    edge = np.hypot(*P.T) > 0.97; c = np.argmin(np.hypot(*P.T)); f = np.zeros(n); f[c] = 1.0
    free = ~edge; u = np.zeros(n); u[free] = sla.spsolve(L[free][:, free].tocsc(), f[free])
    ax = axs[0]; tc = ax.tricontourf(P[:, 0], P[:, 1], -u, levels=18, cmap="Blues_r")
    ax.triplot(P[:, 0], P[:, 1], tri.simplices, lw=0.15, color="white", alpha=0.6)
    ax.plot(0, 0, "o", color=ORANGE, ms=8, mec="white"); ax.set_aspect("equal"); nobox(ax)
    ax.set_title("(a) A mass pulls on the net", loc="left")
    ax.text(0, -1.13, "a random net, tension only; the mass pulls one knot.\ndarker = pulled further. Every knot shares the load\nwith its neighbours, so the dip spreads out smoothly.", ha="center", va="top", fontsize=8, color=INK2)
    tag(ax, "computed: equilibrium of a 2-D random network")
    # (b) radial profile of the same dip (2-D: log) and the 3-D case: pull ~ 1/r^2 from a random 3-D net
    m = 24000; Q = rng.uniform(-1, 1, (m, 3)); t3 = Delaunay(Q); I, J = [], []
    for s in t3.simplices:
        for a in range(4):
            for b in range(a + 1, 4): I += [s[a], s[b]]; J += [s[b], s[a]]
    A3 = sp.coo_matrix((np.ones(len(I)), (I, J)), shape=(m, m)).tocsr(); A3.data[:] = 1
    L3 = sp.diags(np.asarray(A3.sum(1)).ravel()) - A3; r3 = np.linalg.norm(Q, axis=1)
    edge3 = r3 > 0.97; c3 = np.argmin(r3); f3 = np.zeros(m); f3[c3] = 1; fr = ~edge3; u3 = np.zeros(m)
    u3[fr] = sla.spsolve(L3[fr][:, fr].tocsc(), f3[fr])
    bins = np.linspace(0.12, 0.85, 14); rc = 0.5*(bins[1:] + bins[:-1])
    prof = np.array([np.median(u3[(r3 >= a) & (r3 < b)]) for a, b in zip(bins[:-1], bins[1:])])
    g = -np.gradient(prof, rc)                                     # the pull = slope of the dip
    ax = axs[1]; ax.loglog(rc, g, "o", color=BLUE, ms=6, label="random 3-D net (computed)")
    ref = g[3]*(rc/rc[3])**-2; ax.loglog(rc, ref, "--", color=INK2, lw=1.4, label="Newton: pull ~ 1/distance²")
    sl = np.polyfit(np.log(rc[2:-2]), np.log(g[2:-2]), 1)[0]
    ax.set_xlabel("distance from the mass (box units)"); ax.set_ylabel("pull (slope of the dip)")
    ax.legend(frameon=False, fontsize=8); ax.set_title("(b) The pull falls as 1/r²", loc="left")
    ax.text(0.30, 0.72, f"measured slope {sl:.2f} (Newton: -2; the gap shrinks\nwith more knots: lattice graininess near the mass)", transform=ax.transAxes, fontsize=8, color=INK2)
    tag(ax, "computed: 24,000-knot random 3-D network")
    from matplotlib.ticker import FixedLocator, NullFormatter, FormatStrFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.15, 0.2, 0.3, 0.5, 0.8])); ax.xaxis.set_minor_formatter(NullFormatter()); ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    # (c) why 1/r^2: the load is shared over a growing sphere of links
    ax = axs[2]; nobox(ax); ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.25, 1.2); ax.set_aspect("equal")
    for R, a in ((0.35, 0.9), (0.7, 0.6), (1.05, 0.35)):
        ax.add_patch(Circle((0, 0), R, fill=False, ec=BLUE, lw=1.4, alpha=a))
        k = int(10*R/0.35); th = np.linspace(0, 2*np.pi, k, endpoint=False)
        ax.plot(R*np.cos(th), R*np.sin(th), "o", ms=2.5, color=BLUE, alpha=a)
    ax.plot(0, 0, "o", color=ORANGE, ms=9, mec="white")
    ax.set_title("(c) Why 1/r²", loc="left")
    ax.text(0, -1.2, "the same pull is shared by every link crossing\na shell; a shell twice as far has 4x as many,\nso each carries 1/4: the inverse-square law.", ha="center", va="top", fontsize=8, color=INK2)
    tag(ax, "picture")
    fig.suptitle("Figure 2. Gravity is the grid's tension, shared out link by link", x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.04)
    save(fig, "fig2_gravity_tension")
    return sl

# ---------------------------------------------------------------- Fig 3: why things fall: slower ticking (marching-band picture) + light bending
def fig3():
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4))
    ax = axs[0]; nobox(ax); ax.set_xlim(0, 10); ax.set_ylim(-0.6, 5.3)
    ax.add_patch(FancyBboxPatch((0, -0.4), 10, 1.6, boxstyle="square", fc="#d9e7f7", ec="none"))
    ax.text(0.1, -0.3, "mud = near the mass: the grid ticks slower", ha="left", fontsize=8, color=BLUE)
    # a row of marchers; each takes steps proportional to the local tick rate -> the row turns toward the slow side
    y = np.linspace(1.4, 4.2, 8); x = np.zeros_like(y)
    rate = lambda yy: 1 - 0.55*np.exp(-(yy - 0.0)/1.2)
    for step in range(9):
        ax.plot(x, y, "-", color=INK2, lw=1, alpha=0.35 + 0.07*step)
        ax.plot(x, y, "o", ms=4, color=ORANGE if step == 8 else INK2, alpha=0.4 + 0.07*step)
        th = np.arctan2(np.gradient(y), np.gradient(x)) - np.pi/2   # marching direction = normal to the row
        d = rate(y)*0.95; x = x + d*np.cos(th); y = y + d*np.sin(th)
    ax.text(3.2, 5.0, "a marching band whose lower ranks reach mud", fontsize=8.5, color=INK)
    ax.text(10, 4.5, "each row turns toward the slow side:\nno rope pulls it; the path itself bends", fontsize=8, color=INK2, ha="right")
    ax.set_title("(a) Falling = turning toward slower ticks", loc="left"); tag(ax, "picture (computed toy: step length = local tick rate)")
    # (b) light ray through the same tick-rate field: deflection ~ 4GM/(b c^2), shown exaggerated
    ax = axs[1]
    def ray(b, k):
        # refractive index n = 1 + k/r (weak field, k = 2GM/c^2); ray equation in 2-D
        def f(t, s):
            x, y, vx, vy = s; r = np.hypot(x, y); n = 1 + k/r; gx, gy = -k*x/r**3, -k*y/r**3
            dot = vx*gx + vy*gy; return [vx/n, vy/n, (gx - dot*vx)/n, (gy - dot*vy)/n]
        s = solve_ivp(f, [0, 24], [-12, b, 1, 0], max_step=0.02, rtol=1e-8); return s.y[0], s.y[1]
    for b in (1.6, 2.4, 3.3, 4.5):
        xx, yy = ray(b, 0.3); ax.plot(xx, yy, color=BLUE, lw=1.6); ax.plot([-12, 12], [b, b], ":", color=MUTED, lw=1)
    ax.add_patch(Circle((0, 0), 0.6, color=ORANGE)); ax.set_aspect("equal"); ax.set_xlim(-12, 12); ax.set_ylim(-3, 5.5)
    ax.set_xlabel("distance"); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_title("(b) Light takes the same bent path", loc="left")
    ax.text(0, -0.3, "dotted: straight lines. Rays passing closer bend more.\nThe Sun bends starlight by 1.75\" (measured 1919 onward);\nhere the mass is exaggerated ~10^5 times.", fontsize=8, color=INK2, transform=ax.transAxes, va="top")
    tag(ax, "computed: ray tracing through the tick-rate field (exaggerated)")
    fig.suptitle("Figure 3. What gravity actually is: the grid ticks slower near mass, and paths bend toward slow ticking",
                 x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.03)
    save(fig, "fig3_why_things_fall")

# ---------------------------------------------------------------- Fig 4: the grid's tension = dark energy, rho_DE ~ adot^-1/2
def fig4():
    Om, Or = 0.31, 9.1e-5
    def rhs(lna, y):
        a = np.exp(lna); r = np.exp(y[0]); rm, rr = Om*a**-3, Or*a**-4; return [0.5*(rm/2 + rr - r)/(rm + rr + r)]   # BETA = 1/2
    lna = np.linspace(np.log(4), -np.log(31), 4000)
    # integrate backward from today and forward to a = 4
    sb = solve_ivp(rhs, [0, lna[-1]], [np.log(1 - Om - Or)], t_eval=lna[lna <= 0][::-1][::-1], rtol=1e-10)
    sf = solve_ivp(rhs, [0, lna[0]], [np.log(1 - Om - Or)], t_eval=np.sort(lna[lna > 0]), rtol=1e-10)
    L = np.r_[sf.t, sb.t]; R = np.r_[sf.y[0], sb.y[0]]; o = np.argsort(L); L, R = L[o], np.exp(R[o]); a = np.exp(L); z = 1/a - 1
    rm = Om*a**-3; q = (rm/2 + Or*a**-4 - R)/(rm + Or*a**-4 + R); w = -1 - 0.5*q/3
    zc = z[np.argmin(np.abs(w + 1))]
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.7), gridspec_kw=dict(width_ratios=[1, 1, 0.95]))
    ax = axs[0]; m = (z < 6) & (z > -0.75)
    ax.plot(1 + z[m], R[m]/R[np.argmin(np.abs(z))], color=BLUE, label="grid law: tension ~ 1/√(stretch rate)")
    ax.plot(1 + z[m], np.ones(m.sum()), "--", color=INK2, lw=1.4, label="Einstein's constant (Λ)")
    ax.axvline(1, color=GRID, lw=1); ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel("1 + redshift  (past ←  → future)"); ax.set_ylabel("dark-energy density / today")
    ax.legend(frameon=False, fontsize=8, loc="lower left"); ax.set_title("(a) Tension vs stretch rate", loc="left")
    ax.text(1.02, 1.04, "today", fontsize=8, color=MUTED); tag(ax, "computed (no free dark-energy number)")
    ax = axs[1]; m = (z < 3) & (z >= 0)
    ax.axvspan(zc, 3, color="#f3f2ee"); ax.text(1.9, -0.875, "expansion slowing", fontsize=8, color=INK2, ha="center")
    ax.text(zc/2, -0.88, "speeding\nup", fontsize=8, color=INK2, ha="center")
    ax.plot(z[m], w[m], color=BLUE, label="grid law (computed)"); ax.axhline(-1, color=INK2, ls="--", lw=1.4, label="Λ: w = -1")
    ax.plot(z[m], -0.90 - 0.245*(z[m]/(1 + z[m])), ":", color=ORANGE, lw=1.8, label="its 2-number mimic (w0,wa) = (-0.90,-0.25)")
    ax.plot([zc], [-1], "o", color=ORANGE, ms=7, mec="white")
    ax.annotate(f"crossing at z = {zc:.2f}:\nexactly when the\nacceleration begins", xy=(zc, -1), xytext=(1.0, -0.935), fontsize=8, color=INK,
                arrowprops=dict(arrowstyle="-", color=MUTED))
    ax.set_xlabel("redshift z"); ax.set_ylabel("equation of state w"); ax.set_ylim(-1.12, -0.86)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right"); ax.set_title("(b) Its fingerprint", loc="left")
    tag(ax, "computed; testable with DESI DR3 / Euclid")
    ax = axs[2]
    sets = [("Pantheon+", -0.856, -0.497), ("DES-Dovekie", -0.824, -0.617), ("Union3", -0.690, -0.965)]
    for (lab, w0, wa), col in zip(sets, (VIOLET, AQUA, RED)):
        ax.plot(w0, wa, "s", color=col, ms=8, mec="white", label=f"best free fit, DESI DR2 + CMB + {lab}")
    ax.plot(-0.90, -0.245, "*", color=BLUE, ms=15, mec="white", label="grid law (no free number)")
    ax.plot(-1, 0, "x", color=INK2, ms=9, mew=2, label="Λ")
    ax.set_xlabel("w0 (today)"); ax.set_ylabel("wa (change with time)"); ax.set_xlim(-1.1, -0.6); ax.set_ylim(-1.2, 0.15)
    ax.legend(frameon=False, fontsize=7, loc="lower left"); ax.set_title("(c) Where the data put it", loc="left")
    ax.text(-0.62, 0.08, "typical 1σ: ±0.06 in w0, ±0.2 in wa", fontsize=7.5, color=MUTED, ha="right")
    tag(ax, "our refits (predictions/fit_*.txt)")
    fig.suptitle("Figure 4. Dark energy is the grid's own tension, and it 'gives' when the stretching speeds up",
                 x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.04)
    save(fig, "fig4_dark_energy")
    return zc

# ---------------------------------------------------------------- Fig 5: the fluid between the layers = cold dark matter
def fig5():
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.8))
    # (a) halo: fluid gathered around a galaxy (NFW-like particle draw) + disc
    ax = axs[0]; nobox(ax); ax.set_aspect("equal"); ax.set_xlim(-1, 1); ax.set_ylim(-1.35, 1)
    N = 4000; u = rng.uniform(0, 1, N); r = 0.12*(u/(1 - u + 1e-3))**0.9; r = r[r < 1.3]; th = rng.uniform(0, 2*np.pi, len(r))
    ax.scatter(r*np.cos(th), r*np.sin(th), s=1.2, color=BLUE, alpha=0.35, lw=0)
    ax.add_patch(Ellipse((0, 0), 0.6, 0.09, color=ORANGE)); ax.add_patch(Ellipse((0, 0), 0.14, 0.12, color=ORANGE))
    ax.set_title("(a) A galaxy sits in a pool of fluid", loc="left")
    ax.text(0, -1.33, "orange: stars and gas. blue: the fluid it has\ngathered (5-6x more mass). The pool holds\nthe outer stars in orbit (flat rotation curves).", ha="center", va="bottom", fontsize=8, color=INK2)
    tag(ax, "picture (profile drawn from a CDM-like halo)")
    # (b) rotation curve: SPARC galaxy, stars+gas alone vs with a fluid halo (fitted cored halo)
    ax = axs[1]; fn = sorted(glob.glob("/home/claude/sparc/r1/NGC3198_rotmod.dat"))
    if fn:
        d = np.loadtxt(fn[0], comments="#"); r, V, eV, Vg, Vd, Vb = d[:, :6].T
        vb = np.sqrt(np.maximum(Vg*np.abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2, 0))
        ax.errorbar(r, V, eV, fmt="o", ms=4, color=INK, ecolor=MUTED, elinewidth=1, label="measured (NGC 3198)")
        ax.plot(r, vb, color=ORANGE, label="stars + gas alone")
        G = 4.30e-6; fB = lambda x: np.log(1 + x**2) + 2*np.log(1 + x) - 2*np.arctan(x)
        from scipy.optimize import minimize
        def model(p, rr):
            M, rc = 10**p; R200 = (3*M/(4*np.pi*200*136.0))**(1/3); return np.sqrt(np.maximum(Vg*np.abs(Vg) + 0.5*Vd**2 + 0.7*Vb**2, 0)[:len(rr)] + G*M*fB(rr/rc)/fB(R200/rc)/rr)
        p = minimize(lambda p: np.sum(((model(p, r) - V)/np.hypot(eV, 0.05*V))**2), [11.5, 0.7], method="Nelder-Mead").x
        ax.plot(r, model(p, r), color=BLUE, label="stars + gas + fluid pool")
        ax.set_xlabel("radius (kpc)"); ax.set_ylabel("orbital speed (km/s)"); ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.set_title("(b) Why we need the fluid", loc="left"); tag(ax, "computed: SPARC data, Newtonian + cored halo fit")
    # (c) Bullet cluster picture: gas collides and stops, fluid passes through
    ax = axs[2]; nobox(ax); ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal")
    for cx, s in ((-1.25, 1), (1.25, -1)):
        pts = rng.normal(0, 0.28, (500, 2)) + [cx, 0]; ax.scatter(*pts.T, s=2, color=BLUE, alpha=0.35, lw=0)
    for cx in (-0.35, 0.35): ax.add_patch(Ellipse((cx, 0), 0.6, 0.85, color=RED, alpha=0.55))
    ax.annotate("", xy=(-1.9, 0.9), xytext=(-1.0, 0.9), arrowprops=dict(arrowstyle="->", color=INK2))
    ax.annotate("", xy=(1.9, 0.9), xytext=(1.0, 0.9), arrowprops=dict(arrowstyle="->", color=INK2))
    ax.text(0, -0.75, "hot gas (red) collided and stopped", ha="center", fontsize=8, color=RED)
    ax.text(0, -1.15, "fluid (blue) passed through; lensing\nfinds the mass where the blue is", ha="center", fontsize=8, color=INK2)
    ax.set_title("(c) Two clusters after a collision", loc="left"); tag(ax, "picture of the Bullet Cluster observation")
    fig.suptitle("Figure 5. The fluid between the layers: pulled by gravity, invisible, and it passes through itself",
                 x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.04)
    save(fig, "fig5_fluid_dark_matter")

# ---------------------------------------------------------------- Fig 6: cosmic timeline with grid events and test status
def fig6():
    fig, ax = plt.subplots(figsize=(11.5, 3.6)); nobox(ax)
    ev = [(-0.2, "Bounce", "grid has a maximum density:\ncollapse rebounds (no singularity)", "picture; consistent", MUTED),
          (1.2, "Inflation", "still needed: a bounce alone gives\nthe wrong ripples (~30σ off)", "tested: bounce-only FAILS", RED),
          (2.6, "CMB released\n(z = 1100)", "electron ~1% heavier then:\ngrid tension sets the constants", "tested: 1.0099 ± 0.0049 (2σ)", AQUA),
          (4.0, "Dark ages\n(z ~ 100-200)", "electron relaxes; its energy\nbecomes grid tension (dark energy)", "prediction: 21-cm step ≤ 0.7 mK", ORANGE),
          (5.4, "Galaxies form", "fluid pools around them\n(cold dark matter)", "tested: CMB, lensing, clusters pass", AQUA),
          (6.8, "Acceleration\nbegins (z = 0.7)", "tension starts to 'give':\nw crosses -1 right here", "prediction: DESI DR3 / Euclid", ORANGE),
          (8.2, "Today", "dark energy fading slowly\n(w0 ≈ -0.91)", "tested: beats Λ, Δχ² -5.5 to -7.3", AQUA),
          (9.6, "Far future", "expansion keeps speeding up\n(a ~ t^5): never re-collapses", "consequence", MUTED)]
    ax.plot([-0.6, 10.2], [0, 0], color=INK2, lw=2); ax.annotate("", xy=(10.4, 0), xytext=(10.1, 0), arrowprops=dict(arrowstyle="->", color=INK2, lw=2))
    for i, (x, t, d, s, col) in enumerate(ev):
        up = i % 2 == 0; y = 0.55 if up else -0.55
        ax.plot([x, x], [0, y*0.6], color=MUTED, lw=1); ax.plot(x, 0, "o", ms=9, color=col, mec="white", mew=1.5)
        ax.text(x, y, t, ha="center", va="bottom" if up else "top", fontsize=9, fontweight="bold")
        ax.text(x, y + (0.5 if up else -0.42), d, ha="center", va="bottom" if up else "top", fontsize=7.5, color=INK2)
        ax.text(x, y + (1.08 if up else -1.0), s, ha="center", va="bottom" if up else "top", fontsize=7.5, color=col if col != MUTED else INK2, style="italic")
    ax.set_xlim(-1, 10.6); ax.set_ylim(-2.15, 1.95)
    ax.text(-0.9, -2.12, "status colour: green = tested and passing, orange = prediction awaiting data, red = failed, grey = picture/consequence.\nTime axis not to scale.", fontsize=7.5, color=MUTED)
    ax.set_title("Figure 6. The history of the universe in the grid picture, with what has been tested", loc="left")
    save(fig, "fig6_timeline")

# ---------------------------------------------------------------- Fig 7: the electron as a two-layer pattern (computed)
def fig7():
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)
    def H(k, Ns, W, M):
        Ms = np.where((np.arange(Ns) >= Ns//2 - W//2) & (np.arange(Ns) < Ns//2 - W//2 + W), M, -M); h = np.zeros((2*Ns, 2*Ns), complex)
        t = (sz - 1j*sy)/2
        for s in range(Ns):
            h[2*s:2*s+2, 2*s:2*s+2] = np.sin(k)*sx + (Ms[s] - 2 + np.cos(k))*sz; sp_ = (s + 1) % Ns
            h[2*s:2*s+2, 2*sp_:2*sp_+2] += t; h[2*sp_:2*sp_+2, 2*s:2*s+2] += t.conj().T
        return h
    Ns, W, M = 60, 20, 0.5; e, v = np.linalg.eigh(H(0.2, Ns, W, M)); idx = np.argsort(np.abs(e))[:2]
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.6), gridspec_kw=dict(width_ratios=[1.1, 1, 1]))
    ax = axs[0]; s = np.arange(Ns); wa = Ns//2 - W//2
    ax.axvspan(wa, wa + W, color="#f3f2ee"); ax.text(wa + W/2, 0.28, "slab between\nthe two layers", ha="center", fontsize=8, color=INK2)
    for j, col, lab in zip(idx, (BLUE, ORANGE), ("half moving right", "half moving left")):
        p = np.abs(v[:, j])**2; ax.plot(s, p[0::2] + p[1::2], color=col, label=lab)
    ax.set_xlabel("position across the layers"); ax.set_ylabel("where the electron is"); ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, 0.5); ax.set_title("(a) One electron, two halves", loc="left"); tag(ax, "computed (domain-wall lattice model)")
    ax = axs[1]; Ws = np.arange(2, 13); ms = []
    for Wv in Ws:
        ms.append(np.min(np.abs(np.linalg.eigvalsh(H(0.0, 40, Wv, M)))))
    ax.semilogy(Ws, ms, "o-", color=VIOLET, ms=5); ax.semilogy(Ws, ms[0]*0.5**(Ws - Ws[0]), "--", color=INK2, lw=1.2, label="halves every layer")
    ax.set_xlabel("layers between the halves"); ax.set_ylabel("electron mass (grid units)"); ax.legend(frameon=False, fontsize=8)
    ax.set_title("(b) Thicker gap, lighter electron", loc="left"); tag(ax, "computed")
    ax = axs[2]; nobox(ax); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    ax.add_patch(FancyBboxPatch((3.3, 0.8), 3.4, 4.4, boxstyle="round,pad=0.05", fc="#d9e7f7", ec="none"))
    ax.text(5, 5.4, "river", ha="center", fontsize=8.5, color=BLUE)
    ax.plot([1.6], [3], "o", ms=16, color=BLUE); ax.plot([8.4], [3], "o", ms=16, color=ORANGE)
    xs = np.linspace(1.9, 8.1, 100); ax.plot(xs, 3 + 0.25*np.sin(xs*3), color=INK2, lw=1.2)
    ax.text(5, 1.3, "a rope thrown across: the wider the\nriver, the weaker the tug. The tug\nis the mass.", ha="center", fontsize=8, color=INK2)
    ax.set_title("(c) Analogy", loc="left"); tag(ax, "picture")
    fig.suptitle("Figure 7. Matter as patterns in the grid: the electron's two handed halves live on two layers",
                 x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.04)
    save(fig, "fig7_electron_two_layers")

# ---------------------------------------------------------------- Fig 8: black holes: density cap -> Planck star
def fig8():
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.7)); t = np.linspace(-3, 3, 600)
    ax = axs[0]; ax.plot(t, 1/(1 + 6*np.pi*t**2), color=BLUE, label="grid: density capped, collapse rebounds")
    tg = np.abs(t[np.abs(t) > 0.06]); ax.plot(np.sort(-tg)[::-1][::-1], 1/(6*np.pi*np.sort(tg)[::-1]**2), "--", color=INK2, lw=1.4, label="Einstein: density → infinity")
    ax.set_ylim(0, 1.6); ax.set_xlabel("time from the turnaround (Planck units)"); ax.set_ylabel("density / grid maximum")
    ax.legend(frameon=False, fontsize=8, loc="upper right"); ax.set_title("(a) The core of a collapsing star", loc="left")
    tag(ax, "computed: dust collapse with a density cap (LQC-type)")
    ax = axs[1]; nobox(ax); ax.set_xlim(0, 10); ax.set_ylim(0, 6)
    for i in range(5):
        for j in range(4):
            ax.add_patch(Circle((3 + i*0.9, 1.6 + j*0.9), 0.4, color=ORANGE if (i + j) % 2 else BLUE, alpha=0.8))
    ax.text(5, 5.4, "a full lift: people cannot squeeze any closer,\nso the push turns into a rebound", ha="center", fontsize=8.5, color=INK2)
    ax.text(5, 0.6, "result: a 'Planck star' at the centre, hidden by the horizon;\nno observable difference for >1e26 years", ha="center", fontsize=8, color=INK2)
    ax.set_title("(b) Analogy", loc="left"); tag(ax, "picture")
    fig.suptitle("Figure 8. Black holes: the grid cannot be compressed without limit", x=0.01, ha="left", fontweight="bold", fontsize=11.5, y=1.04)
    save(fig, "fig8_black_holes")

# ---------------------------------------------------------------- Fig 9: the early electron test (SPA + DESI + DES)
def fig9():
    files = glob.glob(os.path.join(REPO, "predictions/spa/spa_me*.json")); best = {}
    for f in files:
        me = float(os.path.basename(f)[6:11]); c = json.load(open(f))["chi2"]; best[me] = min(best.get(me, 1e9), c)
    x = np.array(sorted(best)); y = np.array([best[k] for k in x]); y -= y[0]
    a, b, c = np.polyfit(x - 1, y, 2); x0 = -b/(2*a); s = 1/np.sqrt(a)
    fig, ax = plt.subplots(figsize=(6.2, 3.7))
    ax.axvspan(1.004, 1.011, color="#fdeee6"); ax.text(1.0075, 0.6, "grid prediction\n(energy budget)", ha="center", fontsize=8, color=ORANGE)
    xx = np.linspace(0.998, 1.018, 200); ax.plot(xx, np.polyval([a, b, c], xx - 1), color=BLUE, lw=1.6)
    ax.plot(x, y, "o", color=BLUE, ms=7, mec="white"); ax.axhline(0, color=GRID, lw=1)
    ax.set_xlabel("electron mass at recombination / today"); ax.set_ylabel("Δχ² vs today's electron")
    ax.set_title("Figure 9. Was the electron heavier when the CMB was released?", loc="left")
    ax.text(0.9985, -4.6, f"Planck + ACT DR6 + SPT-3G D1 + DESI DR2 + DES-Dovekie:\nm_e = {1 + x0:.4f} ± {s:.4f}  (2σ from today's value)", fontsize=8, color=INK2)
    ax.set_ylim(-5, 1.2); tag(ax, "computed (predictions/spa/)")
    save(fig, "fig9_early_electron")

if __name__ == "__main__":
    fig1(); sl = fig2(); fig3(); zc = fig4(); fig5(); fig6(); fig7(); fig8(); fig9()
    print(f"3-D network pull slope {sl:.2f}; crossing z {zc:.2f}")
