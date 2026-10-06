"""When did dark matter stop trading heat with the hot early universe? Time-temperature windows (radiation era:
t[s] ~ 2.42 / sqrt(g*) / (T/MeV)^2). Limits marked (mem) are from memory of the literature, not looked up here."""
import numpy as np
def t_of_T(T_MeV, g):
    return 2.42/np.sqrt(g)/T_MeV**2
def fmt(t):
    for u, s in ((3.156e7, "yr"), (86400, "days"), (3600, "h"), (1, "s")):
        if t >= u: return f"{t/u:.1f} {s}"
    return f"{t:.0e} s"
rows = []
# thermal relic ('freeze-out'): decouples at T ~ m/20
for m, lab in ((10, "10 MeV"), (1e3, "1 GeV"), (1e5, "100 GeV"), (1e7, "10 TeV")):
    T = m/20; g = 10.75 if T < 100 else 106.75
    rows.append((f"thermal relic, mass {lab}: freezes out at T ~ m/20", T, fmt(t_of_T(T, g))))
# lightest thermal relic allowed by BBN / N_eff (mem: ~10 MeV)
# warm limit: Lyman-alpha thermal-relic mass > 5.7 keV (mem, this project's ledger) -> non-relativistic at T ~ m/3
T = 5.7e-3/3; rows.append(("latest it can become 'cold' (Lyman-alpha, mass > 5.7 keV, mem)", T, fmt(t_of_T(T, 3.36))))
T = 1e-3*1e-3*0.26; rows.append(("matter overtakes radiation (structure starts growing)", 8.0e-7, "~50,000 yr"))
print(f"{'event':72s} {'temperature (MeV)':>18s} {'time after start':>18s}")
for lab, T, t in rows: print(f"{lab:72s} {T:18.2e} {t:>18s}")
print("\nRequired by data (mem): if dark matter was EVER in heat contact (a thermal relic):")
print("  - it must be heavier than ~10 MeV (else it upsets Big Bang element-making / extra radiation) -> contact ended before ~3 s")
print("  - its amount (5.4x ordinary matter) fixes its annihilation strength <sigma v> ~ 3e-26 cm^3/s ('WIMP' case), any mass 10 MeV-100 TeV")
print("  - heat contact with ordinary matter/light must be tiny afterwards (CMB limits on dark-matter scattering)")
print("If it was NEVER in heat contact (axion-like field, primordial black holes, or a dust tied to the grid's clock): no 'when' exists;")
print("  its coldness is set at birth.")
