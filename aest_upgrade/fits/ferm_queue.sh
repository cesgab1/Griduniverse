cd /home/claude/fullfit
for W in 0 0.002 0.005; do OMP_NUM_THREADS=1 EVAL=1 AEST_BETA=0.5 AEST_EX=yes FERM_OMEGA=$W python3 run_pact_class.py ferm_eval_$W '{"H0": 68.502588, "ns": 0.97434898, "ombh2": 0.022548605, "omch2": 0.11788114, "tau": 0.059290191, "A_act": 1.0004043, "logA": 3.0495784, "P_act": 1.0}' > pact_class/log_ferm_$W.txt 2>&1; done &
for W in 0.01 0.02 0.05; do OMP_NUM_THREADS=1 EVAL=1 AEST_BETA=0.5 AEST_EX=yes FERM_OMEGA=$W python3 run_pact_class.py ferm_eval_$W '{"H0": 68.502588, "ns": 0.97434898, "ombh2": 0.022548605, "omch2": 0.11788114, "tau": 0.059290191, "A_act": 1.0004043, "logA": 3.0495784, "P_act": 1.0}' > pact_class/log_ferm_$W.txt 2>&1; done &
wait
