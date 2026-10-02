#!/bin/bash
cd /home/claude/fullfit
run(){ me=$1; src=$2; tag=$3; st=$(python3 -c "import json;print(json.dumps(json.load(open('spa/$src.json'))['params']))"); VC_ME=$me python3 run_spa.py $tag "$st" > spa/log_$tag.txt 2>&1; }
(run 1.008 spa_me1.012 spa_me1.008b; run 1.000 spa_me1.004 spa_me1.000b; run 1.016 spa_me1.012 spa_me1.016b) &
(run 1.008 spa_me1.004 spa_me1.008c; run 1.004 spa_me1.012 spa_me1.004b; run 1.012 spa_me1.016 spa_me1.012b) &
wait
