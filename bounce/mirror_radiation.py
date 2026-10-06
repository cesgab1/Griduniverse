"""Supplement (Oct 5 2026): in the mirror-image contraction, modes that entered OUR horizon in the radiation era (k > k_eq)
left the horizon during a RADIATION-dominated contraction. Wands duality: n_s - 1 = 12 w / (1 + 3w) for a contracting fluid
-> w = 1/3 gives n_s = 3 (strongly blue). Also: how far the bounce's own scale lies below what we can observe."""
import numpy as np
k_eq = 0.0104            # 1/Mpc, matter-radiation equality
k_max = 0.2              # 1/Mpc, smallest scales in the CMB fits
w = 1/3; ns_rad = 1 + 12*w/(1 + 3*w)
pred = (k_max/k_eq)**(ns_rad - 1); obs = (k_max/k_eq)**(0.9649 - 1)
print(f"radiation-contraction tilt n_s = {ns_rad:.1f}; power at k = {k_max} vs k_eq: predicted x{pred:.0f}, measured x{obs:.2f}")
lP = 1.616e-35; patch = 0.017e-3                      # visible universe at the cap (reverse_bang/)
print(f"bounce curvature scale ~ Planck length; visible universe at the bounce {patch:.1e} m -> bounce features sit at scales "
      f"{patch/lP:.0e} times smaller than the largest we can observe (unless inflation stretched them)")
