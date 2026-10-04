"""ITERATION 70 (pre-registered in PREREG_70.md)."""
import numpy as np
base = 17.439          # iteration 68: plain Planck length, causal diamond (past light cone), factor 1/2 -> 17.44 x measured
out = ["ITERATION 70: which smallest length? (expectations pre-registered; P_E and P_0 not blind)", "",
       " principle                         l^2 / l_P^2   dark energy / measured   entropy per element (black holes)"]
for name, f in (("P_0 plain Planck", 1.0), ("P_S black-hole entropy (Hawking)", 4.0), ("P_E Einstein coupling (8 pi)", 8*np.pi), ("P_A action (16 pi)", 16*np.pi)):
    out.append(f" {name:34s} {f:8.2f}        {base/f:8.2f}                {f/4:6.2f} nats")
out += ["", "Only P_S gives black holes exactly 1 nat per element; it gives dark energy 4.4 x measured.",
        "P_E gives dark energy 0.69 x but needs 2 pi nats per element for black holes. The two natural requirements disagree:",
        "GR + hbar do not single out one smallest length -> dark energy's size from this route: 0.35 - 4.4 x measured (not pinned)."]
txt = "\n".join(out); print(txt); open("iter70_smallest_length.txt", "w").write(txt + "\n")
