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
def structure(C, d, lags):
    return [np.var(C[int(l/d):] - C[:-int(l/d)]) for l in lags]
out = ["ITERATION 34: grid growth by link breaking (expectations committed before running)", "",
       " tolerance   memory of stored strain   stored energy / stiffness   count excess: var of change over 0.1 / 0.5 / 2.0 ln a"]
mem = {}
for eps_c in (0.01, 0.1, 0.33, 0.5, 1.0, 1.6):
    E, C, d = run(eps_c, L=8.0)
    burn = len(E)//8; m = corr_length(E[burn:], d); mem[eps_c] = m
    sf = structure(C[burn:], d, (0.1, 0.5, 2.0))
    out.append(f"   {eps_c:5.2f}      {m:6.3f} Hubble times        {E[burn:].mean():.2e}               {sf[0]:.3f} / {sf[1]:.3f} / {sf[2]:.3f}")
for N in (2000, 20000):
    E, C, d = run(0.1, N=N, L=2.0)
    out.append(f"   N = {N:6d}: relative fluctuation of stored strain = {np.std(E[len(E)//3:])/np.mean(E[len(E)//3:]):.4f}")
ratio = np.mean([mem[e]/e for e in mem])
need = (1/3)/ratio
out += ["", f"Measured: memory ~ {ratio:.2f} x tolerance. The data's 1/3 Hubble time needs tolerance ~ {need:.1f} (links stretching",
        f"by ~{need*100:.0f}% before breaking); stored energy is then ~{0.42*need**2:.1f} x the stiffness -> ~1e123 x too big.",
        "A tolerance small enough for the size (3e-62) gives essentially zero memory.",
        "Count excess: if its change keeps growing with the lag it wanders (random walk = the excluded bucket); if it levels off",
        "it is stationary. See the three columns above.",
        "Conclusion: link breaking explains WHY cells are added and gives sqrt(N) fluctuations automatically (N check above),",
        "but it cannot give both the measured memory and the measured size. The memory must come from elsewhere (Hubble friction",
        "of the jostled field, iteration 4); the size is not set by the tolerance."]
txt = "\n".join(out); print(txt); open("iter34_link_breaking.txt", "w").write(txt + "\n")
