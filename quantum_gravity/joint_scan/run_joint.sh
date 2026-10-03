#!/bin/bash
# Joint scan of counting exponent BETA and link-stretch exponent s: d ln rho_DE/d ln a = BETA (q + 1 - s); simplified q.
cd /home/claude/griduniverse/cosmology_fits
OUT=/home/claude/griduniverse/quantum_gravity/joint_scan
job() { SN=$1; B=$2; S=$3
  if [ "$B" = "LCDM" ]; then r=$(SNSET=$SN ONLY=LCDM python3 fit_law.py | grep ONLYRESULT); echo "$SN LCDM x $r"
  else r=$(SNSET=$SN BETA=$B STRETCH=$S ONLY=LAW python3 fit_law.py | grep ONLYRESULT); echo "$SN $B $S $r"; fi; }
export -f job
{ for SN in PANTHEON DESY5 UNION3; do echo "$SN LCDM x"; for B in 0.25 0.375 0.5 0.625 0.75 1.0; do for S in 0.25 0.5 0.75 1 1.25 1.5 1.75; do echo "$SN $B $S"; done; done; done; } \
  | xargs -P 2 -L 1 bash -c 'job "$@"' _ > $OUT/joint_raw.txt
touch $OUT/DONE
