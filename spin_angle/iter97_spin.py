import numpy as np
c = 2.998e5; N, Mg, sig, R, vr, Rd = 1000, 1e11, 1000.0, 1000.0, 200.0, 5.0   # km/s, kpc
J_orb = N*Mg*sig*R                 # rough scale of orbital angular-momentum budget (|r x v| summed, no cancellation)
J_spin_aligned = N*Mg*vr*2*Rd      # exponential disc: J = 2 M v R_d
J_spin_random = np.sqrt(N)*Mg*vr*2*Rd
out = [f"spin/orbital angular momentum: aligned {J_spin_aligned/J_orb:.1e}, random {J_spin_random/J_orb:.1e}",
       f"spin energy as extra mass (v_rot/c)^2 = {(vr/c)**2:.1e};  frame-dragging ~ (sigma/c)^2 = {(sig/c)**2:.1e}",
       f"needed extra: x5 in mass = +400% -> spin effects short by a factor {4/((sig/c)**2):.0e} or more",
       "viewing angle: anisotropy beta = 0.5 vs 0 changes the virial mass estimate by up to ~x1.3-2 for one cluster; averages to ~x1"]
print("\n".join(out)); open("iter97_spin.txt", "w").write("\n".join(out) + "\n")
