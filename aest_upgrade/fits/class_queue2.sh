cd /home/claude/fullfit
S=$(python3 -c "import json; d=json.load(open('pact_class/aest_b05_ex.json'))['params']; d2={'H0':d['H0'],'ns':d['n_s'],'ombh2':d['omega_b'],'omch2':d['omch2'],'tau':d['tau_reio'],'A_act':d['A_act'],'logA':d['logA'],'P_act':1.0}; print(json.dumps(d2))")
(OMP_NUM_THREADS=1 AEST_BETA=0.5 AEST_EX=yes AEST_ME=1.0043 python3 run_pact_class.py aest_b05_ex_me043 "$S" > pact_class/log_me043.txt 2>&1) &
(OMP_NUM_THREADS=1 AEST_BETA=0.5 AEST_EX=yes AEST_ME=1.01 python3 run_pact_class.py aest_b05_ex_me100 "$S" > pact_class/log_me100.txt 2>&1) &
wait
