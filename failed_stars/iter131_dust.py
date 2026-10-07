"""PREREG 131: where the dust came from (budget), how it was gathered (snowplough shell), why it stays (drift)."""
import numpy as np
out=[]; P=lambda s="": (print(s), out.append(s))
P("(a) Milky Way dust budget (published inputs, Draine 2009 / Zhukovska 2008, 2016 / Slavin 2015):")
dust_obs = 2.5e7; inj_star = 0.005
for lab, tau in (("Draine 2009 shock lifetime", 4e8), ("Zhukovska 2016 best fit", 3.5e8), ("Slavin 2015 silicates", 2.5e9)):
    made = inj_star*tau
    P(f"   destruction time {tau/1e9:.2f} Gyr ({lab}): stardust standing stock = {made:.1e} Msun = {made/dust_obs*100:.0f}% of "
      f"observed -> {100-min(100,made/dust_obs*100):.0f}% must regrow in clouds")
sn_rate = 0.02; sn_dust = 0.5       # ~2 per century (assumed, memory); SN 1987A 0.5-0.8 Msun (ALMA/Herschel)
P(f"   if every supernova kept SN-1987A-like dust ({sn_dust} Msun, rate ~{sn_rate}/yr, assumed): +{sn_rate*sn_dust:.3f} Msun/yr "
  f"(x{sn_rate*sn_dust/inj_star:.0f} the other stellar sources) -> survival through the reverse shock decides it")
P("\n(b) Gathering: Per-Tau shell snowplough (radius 78 pc, Bialy 2021); ambient density not published -> range")
pc=3.086e18; mH=1.6735e-24; Msun=1.989e33
Vol = 4/3*np.pi*(78*pc)**3
clouds = 1.9e4 + 1.5e4   # Perseus (extinction map) + Taurus (H2)
for n in (1, 2, 5):
    m = Vol*n*1.4*mH/Msun
    P(f"   n = {n} cm^-3: swept gas {m:.1e} Msun (dust ~{m/100:.0f}-{m/200:.0f} Msun); Perseus+Taurus clouds = {clouds/m*100:.0f}% of it")
P("\n(c) Why mm dust stays in disks: classical drag drift (Epstein), disk = 1% of central mass, Sigma ~ 1/r, T ~ r^-1/2")
G=6.674e-8; k=1.381e-16; AU=1.496e13; yr=3.156e7
def drift_time(Mstar, L, Mdisk_frac, rc_AU, a_cm=0.1, rho_s=1.6):
    M = Mstar*Msun; Md = Mdisk_frac*M
    r = 0.5*rc_AU*AU
    Sigc = Md/(2*np.pi*(rc_AU*AU)**2)           # Sigma(r) = Sigc * rc/r inside rc
    Sig = Sigc*(rc_AU*AU)/r
    T = 280*L**0.25*(r/AU)**-0.5
    cs2 = k*T/(2.33*mH); vK = np.sqrt(G*M/r)
    eta = 0.5*cs2/vK**2*2.75
    St = np.pi/2*a_cm*rho_s/Sig
    vr = 2*St/(1+St**2)*eta*vK
    return r/vr/yr, St, eta
for lab, M, L, rc in (("Sun-like star (1 Msun, 1 Lsun, disk 100 AU)",1.0,1.0,100),
                      ("brown dwarf (0.05 Msun, 0.003 Lsun, disk 20 AU)",0.05,0.003,20),
                      ("IC 348-type (0.002 Msun ~2 MJup, 3e-4 Lsun, disk 10 AU)",0.002,3e-4,10)):
    t, St, eta = drift_time(M, L, 0.01, rc)
    P(f"   {lab:58s}: 1-mm grain at half disk radius drifts in {t:.1e} yr (Stokes {St:.2f}, eta {eta:.3f})")
P("   observed: mm grains present in brown-dwarf disks at ~1-3 Myr ages (Ricci 2014; Testi 2022 detections above)")
open("RESULT_131_numbers.txt","w").write("\n".join(out)+"\n")
