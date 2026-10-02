#!/bin/bash
cd /home/claude/fullfit
run(){ me=$1; st=$(python3 -c "import json;p=json.load(open('pact/vc_me1.01.json'))['params'];print(json.dumps(p))"); VC_ME=$me python3 run_spa.py spa_me$me "$st" > spa/log_me$me.txt 2>&1; }
(run 1.000; run 1.008; run 1.016) &
(run 1.004; run 1.012) &
wait
