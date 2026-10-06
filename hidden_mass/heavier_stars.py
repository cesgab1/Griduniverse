"""Coalesce: masses are counted from light (stars) and X-rays (gas); objects could be heavier than their light suggests (hidden dense
cores). How much heavier would the STARS have to be to close the cluster gap, and is that allowed?"""
cases = {"Coma (whole cluster)": dict(gas=2.5e14, stars=0.5e14, needed=1.6e15),
         "Bullet, main galaxies (100 kpc)": dict(gas=5.5e12, stars=0.54e12, needed=2.7e13 + 5.5e12 + 0.54e12),
         "Bullet, sub galaxies (100 kpc)": dict(gas=2.7e12, stars=0.58e12, needed=1.5e13)}
for k, d in cases.items():
    f = (d["needed"] - d["gas"])/d["stars"]
    print(f"{k:32s}: stars would have to weigh x{f:.0f} what their light suggests")
print("""
Independent checks on how heavy stars really are (mem):
  * gas mass is NOT from light: X-ray brightness gives it directly (and the SZ effect agrees)
  * star masses checked by motions in galaxy centres and by strong lensing of elliptical galaxies' inner parts:
    true star mass = 1-2x the light-based estimate (the 'bottom-heavy' debate), not x10-x50
  * star clusters / binaries: individual stellar masses known to ~percent
  * deuterium cap: making stars x27-x50 heavier puts ordinary matter at ~5-6x the Big-Bang value
A star with a hidden neutron star (Thorne-Zytkow) is ~+10% heavier than its light suggests -- far short of x27.""")
