import numpy as np, sys, json
sys.argv = ["x", "0"]
src_code = open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/crack_grid.py").read().split("out = {}")[0]
exec(src_code)
import heapq
def arrival_from(P, adj, src):
    D = {src: 0.0}; h = [(0.0, src)]
    while h:
        d, u = heapq.heappop(h)
        if d > D.get(u, 1e9): continue
        for v, w in adj.get(u, []):
            if d + w < D.get(v, 1e9): D[v] = d + w; heapq.heappush(h, (d + w, v))
    return D
bins = np.radians(np.arange(0, 91, 10))
res = {}
for lab, build in (("(a) horizontal/vertical cracks", lambda: crack_mosaic(3000, True)), ("(b) random-angle cracks", lambda: crack_mosaic(3000, False)),
                   ("(c) random-point grid", lambda: voronoi_grid(3000))):
    prof = [[] for _ in range(len(bins)-1)]
    for rep in range(4):
        P, edges, _ = build(); adj = {}
        for u, v in edges:
            w = np.linalg.norm(P[u]-P[v]); adj.setdefault(u, []).append((v, w)); adj.setdefault(v, []).append((u, w))
        nodes = np.array(list(adj.keys())); cand = nodes[np.linalg.norm(P[nodes]-0.5, axis=1) < 0.15]
        for src in rng.choice(cand, min(12, len(cand)), replace=False):
            D = arrival_from(P, adj, src); ids = np.array(list(D)); t = np.array([D[i] for i in ids]); X = P[ids]-P[src]; r = np.linalg.norm(X, axis=1)
            m = (r > 0.2) & (r < 0.35) & (np.abs(P[ids]-0.5).max(axis=1) < 0.47)
            a = np.abs(np.arctan2(X[m, 1], X[m, 0])); a = np.where(a > np.pi/2, np.pi - a, a)
            q = t[m]/r[m]
            for k in range(len(bins)-1):
                s = (a >= bins[k]) & (a < bins[k+1])
                if s.sum(): prof[k].append(np.median(q[s]))
    med = np.array([np.mean(p) for p in prof]); err = np.array([np.std(p)/np.sqrt(len(p)) for p in prof])
    res[lab] = (med.tolist(), err.tolist())
    axis = (med[0] + med[-1])/2; diag = (med[4])
    print(f"{lab:32s} travel-path/straight-line by direction (0..90 deg in 10-deg bins): " + " ".join(f"{x:.3f}" for x in med))
    print(f"{'':32s} diagonal vs along-axis: {(diag/axis-1)*100:+.1f}%  (+-{np.hypot(err[4], err[0])/axis*100:.1f}%)")
json.dump(res, open("/tmp/claude-0/-home-claude/631bebbd-9799-5238-817c-c2760930ea03/scratchpad/crack_grid_avg.json", "w"))
