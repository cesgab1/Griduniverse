#!/bin/bash
# Scan the link-stretch exponent s (physical link length ∝ a^s) in the law d ln rho_DE/d ln a = (q + 1 - s)/2.
# s = 1: links stretch with space (Claim 1, rho ∝ adot^-1/2).  s = 0: cells added, link length fixed (rho ∝ H^-1/2).
cd /home/claude/griduniverse/cosmology_fits
OUT=/home/claude/griduniverse/quantum_gravity/stretch_scan
for SN in PANTHEON DESY5 UNION3; do
  for Q in simple exact; do
    EX=""; [ $Q = exact ] && EX=1
    SNSET=$SN EXACTQ=$EX ONLY=LCDM python3 fit_law.py | grep ONLYRESULT | sed "s/^/$SN $Q LCDM /" >> $OUT/scan_raw.txt
    for S in -0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.5; do
      SNSET=$SN EXACTQ=$EX STRETCH=$S ONLY=LAW python3 fit_law.py | grep ONLYRESULT | sed "s/^/$SN $Q s=$S /" >> $OUT/scan_raw.txt
    done
  done
done
