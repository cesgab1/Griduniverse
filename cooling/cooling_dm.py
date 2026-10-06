"""Could dark matter be built from cooling (radiation energy lost to expansion)? Energy bookkeeping per comoving volume."""
Om_r, Om_dm = 9.1e-5, 0.26
E_dm_over_Er0 = Om_dm/Om_r
for z in (1090, 3400):
    lost = z            # radiation energy per comoving volume ~ 1/a: lost since z = E_r0 * z
    print(f"radiation energy lost since z = {z}: {lost:.0f} E_r0 = {lost/E_dm_over_Er0:.0%} of dark matter's energy today")
print(f"dark matter / radiation at z = 1090: {Om_dm/(Om_r*1091):.2f} (already ~3x radiation when the glow was released)")
print("temperature contrast at z = 1090: ~1e-5 everywhere -> no 'cooled regions' yet, but dark matter already 5.4x ordinary matter")
