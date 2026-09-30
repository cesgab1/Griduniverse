"""
Grain check for cracked grids (T-junction mosaics made by repeatedly splitting cells).
 (a) cracks only horizontal/vertical (the user's picture)   (b) cracks at random angles   (c) random-point (Voronoi) grid
A signal leaves the centre and travels along the links at a fixed speed. Grain = the arrival time depends on direction.
Measured: path-along-links / straight-line distance vs direction, for points 0.25-0.42 from the centre.
"""
import numpy as np, heapq, json, sys
from scipy.spatial import Voronoi
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def crack_mosaic(ncut, axis_only):
    P = [np.array([0., 0.]), np.array([1., 0.]), np.array([1., 1.]), np.array([0., 1.])]   # node coordinates
    polys = [[0, 1, 2, 3]]; E2P = {}
    key = lambda u, v: (min(u, v), max(u, v))
    for i in range(4): E2P.setdefault(key(i, (i+1) % 4), set()).add(0)
    def perim(pl): return sum(np.linalg.norm(P[pl[i]] - P[pl[(i+1) % len(pl)]]) for i in range(len(pl)))
    per = [perim(polys[0])]
    for _ in range(ncut):
        pi = rng.choice(len(polys), p=np.array(per)/sum(per)); pl = polys[pi]
        th = rng.choice([0, np.pi/2]) if axis_only else rng.uniform(0, np.pi)
        nrm = np.array([np.cos(th), np.sin(th)]); pr = np.array([P[v] @ nrm for v in pl])
        off = rng.uniform(pr.min(), pr.max())
        # find the two crossing edges
        cr = []
        for i in range(len(pl)):
            u, v = pl[i], pl[(i+1) % len(pl)]; du, dv = P[u] @ nrm - off, P[v] @ nrm - off
            if du*dv < 0: cr.append((i, u, v, du/(du - dv)))
        if len(cr) != 2: continue
        newids = []
        for (i, u, v, t) in cr:
            w = len(P); P.append(P[u] + t*(P[v] - P[u])); newids.append(w)
            for q in list(E2P.pop(key(u, v), set())):          # insert w into every cell sharing edge u-v
                ql = polys[q]
                for j in range(len(ql)):
                    if {ql[j], ql[(j+1) % len(ql)]} == {u, v}:
                        ql.insert(j+1, w); break
                E2P.setdefault(key(u, w), set()).add(q); E2P.setdefault(key(w, v), set()).add(q)
        pl = polys[pi]; a, b = newids
        ia, ib = pl.index(a), pl.index(b)
        if ia > ib: ia, ib = ib, ia
        p1 = pl[ia:ib+1]; p2 = pl[ib:] + pl[:ia+1]
        e1 = {key(p1[j], p1[(j+1) % len(p1)]) for j in range(len(p1))}
        e2 = {key(p2[j], p2[(j+1) % len(p2)]) for j in range(len(p2))}
        nq = len(polys); polys[pi] = p1; polys.append(p2); per.append(0)
        for k_ in e2:
            s_ = E2P.setdefault(k_, set())
            if k_ not in e1: s_.discard(pi)
            s_.add(nq)
        for k_ in e1: E2P.setdefault(k_, set()).add(pi)
        per[pi] = perim(p1); per[nq] = perim(p2)
    edges = [k for k, s in E2P.items() if s]
    return np.array(P), edges, len(polys)
def voronoi_grid(npts):
    pts = rng.uniform(-0.2, 1.2, (int(npts*1.96), 2)); vor = Voronoi(pts); V = vor.vertices
    edges = [tuple(r) for r in vor.ridge_vertices if -1 not in r]
    return V, edges, npts
def arrival(P, edges):
    adj = {}
    for u, v in edges:
        w = np.linalg.norm(P[u] - P[v]); adj.setdefault(u, []).append((v, w)); adj.setdefault(v, []).append((u, w))
    src = int(np.argmin(np.linalg.norm(P - 0.5, axis=1))); D = {src: 0.0}; h = [(0.0, src)]
    while h:
        d, u = heapq.heappop(h)
        if d > D.get(u, 1e9): continue
        for v, w in adj.get(u, []):
            if d + w < D.get(v, 1e9): D[v] = d + w; heapq.heappush(h, (d + w, v))
    return src, D
out = {}
for lab, build in (("(a) horizontal/vertical cracks", lambda: crack_mosaic(1500, True)), ("(b) random-angle cracks", lambda: crack_mosaic(1500, False)),
                   ("(c) random-point grid", lambda: voronoi_grid(1500))):
    P, edges, ncell = build(); src, D = arrival(P, edges)
    ids = np.array(list(D.keys())); t = np.array([D[i] for i in ids]); X = P[ids] - P[src]; r = np.linalg.norm(X, axis=1)
    m = (r > 0.25) & (r < 0.42) & (np.abs(P[ids] - 0.5).max(axis=1) < 0.45)
    ang = np.mod(np.arctan2(X[m, 1], X[m, 0]), np.pi/2); ratio = t[m]/r[m]
    near_axis = ratio[(ang < np.radians(8)) | (ang > np.radians(82))]; diag = ratio[np.abs(ang - np.pi/4) < np.radians(8)]
    c4 = np.mean(ratio*np.cos(4*ang))/np.mean(ratio)*2
    print(f"{lab:32s} cells {ncell:5d} | travel-path/straight-line: along axes {np.median(near_axis):.3f}, diagonal {np.median(diag):.3f}"
          f" -> direction difference {(np.median(diag)/np.median(near_axis)-1)*100:+.1f}% | 4-fold grain amplitude {c4*100:+.1f}%")
    out[lab] = dict(P=P.tolist(), edges=[list(map(int, e)) for e in edges], src=src, D={int(k): float(v) for k, v in D.items()},
                    axis=float(np.median(near_axis)), diag=float(np.median(diag)))
json.dump(out, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/crack_grid.json", "w"))
