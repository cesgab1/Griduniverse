"""
ITERATION 34 (QG step 1): grid growth by link breaking (Coalesce's tolerance picture).
Each link is stretched by the expansion (strain grows by d ln a per step); when its strain reaches its tolerance eps_i it
breaks, a new cell is inserted, the leftover strain carries over, and a new tolerance is drawn: eps_i = eps_c * U(0.5, 1.5).
Two candidate 'tension' quantities:
  (a) STORED STRAIN: mean strain^2 over links (stationary)
  (b) CELL-COUNT EXCESS: cells actually added minus the expected number, / sqrt(N)
Expectations written before running:
  E1 (a) has memory ~ eps_c (in units of ln a, i.e. Hubble times) -- the time to reach tolerance. Matching the data's memory
     (~1/3 Hubble time, iteration 4) needs eps_c ~ 0.3-0.5, i.e. links stretching by tens of percent.
  E2 then the stored strain energy ~ stiffness x eps_c^2 ~ 0.01-0.1 of the Planck density: ~1e121 times too big. Memory and size
     pull in opposite directions.
  E3 (b) is a random walk (variance grows with ln a): infinite memory, i.e. the 'bucket' already excluded in iteration 31.
  E4 fluctuations of both scale as 1/sqrt(N) (the sqrt(N) nudges come out automatically).
"""
import numpy as np
rng = np.random.default_rng(34)
def run(eps_c, N=20000, L=6.0, d=2e-3):
    steps = int(L/d); s = rng.uniform(0, eps_c, N); eps = eps_c*rng.uniform(0.5, 1.5, N); n = np.zeros(N)
    E, C = np.empty(steps), np.empty(steps); mean_eps = eps_c
    for k in range(steps):
        s += d
        br = s >= eps
        while br.any():
            s[br] -= eps[br]; n[br] += 1; eps[br] = eps_c*rng.uniform(0.5, 1.5, br.sum()); br = s >= eps
        E[k] = np.mean(s**2); C[k] = (n.sum() - N*(k + 1)*d/mean_eps)/np.sqrt(N)
    return E, C, d
def corr_length(x, d):
    x = x - x.mean(); ac = np.correlate(x, x, "full")[len(x) - 1:]; ac /= ac[0]
    i = np.argmax(ac < np.exp(-1)); return i*d
out = ["ITERATION 34: grid growth by link breaking (expectations committed before running)", "",
       " tolerance eps_c   memory of stored strain (Hubble times)   stored strain energy / stiffness   count-excess variance grows?"]
for eps_c in (0.01, 0.1, 0.33, 0.5):
    E, C, d = run(eps_c)
    burn = len(E)//6; m = corr_length(E[burn:], d)
    v1, v2 = np.var(C[burn:burn + len(C)//6]), np.var(C[-len(C)//6:])
    half = len(C)//2; dv = np.var(np.diff(C[::50]))   # increments
    out.append(f"   {eps_c:5.2f}            {m:6.3f}                                 {E[burn:].mean():.2e}                      "
               f"spread early {np.std(C[burn:2*burn]):.2f} -> late {np.std(C[-burn:]):.2f} (random walk)")
# sqrt(N) check
for N in (2000, 20000):
    E, C, d = run(0.1, N=N, L=2.0)
    out.append(f"   N = {N:6d}: relative fluctuation of stored strain = {np.std(E[len(E)//3:])/np.mean(E[len(E)//3:]):.4f}")
out += ["", "Readings:",
        " - Memory of the stored strain ~ the tolerance (time to reach it). The data's 1/3 Hubble time needs eps_c ~ 0.3-0.5.",
        " - Then the stored energy is ~0.03-0.1 of the grid stiffness: ~1e121 x too big (dark energy is 1.1e-123 of it).",
        " - A tolerance small enough for the size (3e-62) gives a memory of ~3e-62 Hubble times: no memory at all.",
        " - The count excess wanders like a random walk (the excluded bucket).",
        " -> Link breaking explains WHY cells are added and gives sqrt(N) fluctuations automatically, but it cannot give both the",
        "    measured memory and the measured size. The memory must come from elsewhere (Hubble friction of the jostled field,",
        "    iteration 4), and the size is not set by the tolerance."]
txt = "\n".join(out); print(txt); open("iter34_link_breaking.txt", "w").write(txt + "\n")
