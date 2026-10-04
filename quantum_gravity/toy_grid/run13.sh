#!/bin/bash
cd /home/claude/griduniverse/quantum_gravity/toy_grid
for SN in PANTHEON DESY5 UNION3; do for F in 1.25 1.5 2 3; do echo "$SN $F"; done; done | xargs -P 2 -L 1 bash -c 'SNSET=$0 FRAC=$1 python3 iter13_baseline_fraction.py 2>/dev/null | grep RESULT' > iter13_raw_f_gt1.txt
touch DONE13b
