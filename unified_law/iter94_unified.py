"""Iteration 94: unified galaxy-cluster law (PREREG_94.md)."""
G, MSUN, A0 = 6.674e-11, 1.989e30, 1.2e-10
v4 = lambda M: (G * M * MSUN * A0) ** 0.25 / 1e3                       # km/s
kappa = lambda sig, M: (sig * 1e3) ** 4 / (G * M * MSUN * A0)
L = []
# C3 + C1 on NGC 3198
vo = v4(3e10); k = (10 / 150) ** 2; vr = k ** 0.25 * vo
L += [f"C3 NGC 3198: v_ord = (G M a0)^1/4 = {vo:.1f} km/s (observed ~150)",
      f"C1 NGC 3198: k = (10/150)^2 = {k:.4f} -> formula gives v_rand = k^1/4 v_ord = {vr:.1f} km/s, but the input random motion was 10 km/s",
      f"   v_tot = {(vo**2 + vr**2)**0.5:.1f} km/s.  Self-consistency needs k^2 = k (k = 0 or 1). Coma: v_ord = 0 -> k infinite (undefined).", ""]
# C2 one universal kappa
L.append("C2 one kappa for all random-motion systems: kappa = sigma^4 / (G M_b a0)   [MOND isothermal theory: 4/81 = 0.049]")
for name, s, Ms in [("Fornax dSph", 11.7, [2e7, 3e7, 4e7]), ("Coma", 1000, [3e14, 1.5e14, 2e13])]:
    for M in Ms:
        L.append(f"  {name:11s} sigma {s:6.1f} km/s  M_b {M:.1e} Msun  -> kappa = {kappa(s, M):.3f}")
kF, kC = kappa(11.7, 3e7), kappa(1000, 3e14)
L += ["", f"Coma needs kappa {kC/kF:.1f}x Fornax's (central values).",
      f"Equivalently: with Fornax's kappa, Coma's mass from its speeds = {1e24/(kF*G*MSUN*A0):.2e} Msun vs visible 3e14 -> gap x{kC/kF:.1f}"]
txt = "\n".join(L); print(txt); open("iter94_unified.txt", "w").write(txt + "\n")
