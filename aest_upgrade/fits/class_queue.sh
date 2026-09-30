cd /home/claude/fullfit
S=$(cat pact_class/start.json)
(OMP_NUM_THREADS=1 AEST_BETA=0 python3 run_pact_class.py aest_lambda "$S" > pact_class/log_aest_lambda.txt 2>&1;
 OMP_NUM_THREADS=1 AEST_BETA=0.5 AEST_EX=yes python3 run_pact_class.py aest_b05_ex "$S" > pact_class/log_aest_b05_ex.txt 2>&1) &
(OMP_NUM_THREADS=1 AEST_BETA=0.5 AEST_EX=no python3 run_pact_class.py aest_b05 "$S" > pact_class/log_aest_b05.txt 2>&1) &
wait
