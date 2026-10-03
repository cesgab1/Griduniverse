#!/bin/bash
cd /home/claude/fullfit
run(){ al=$1; st=$(python3 -c "import json;print(json.dumps(json.load(open('spa/spa_me1.000.json'))['params']))"); VC_ALPHA=$al VC_ME=1.0 python3 run_spa.py spa_al$al "$st" > spa/log_al$al.txt 2>&1; }
(run 0.996; run 1.004) &
(run 1.008; run 0.992) &
wait
