#!/bin/bash
# waits for certify_joint.py run 8 to finish, then runs the follow-up and the exact verification
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd C:/Projects/Claude4/lb_joint
until grep -q "^wrote" run_V8.log; do sleep 20; done
../.venv/Scripts/python analyze_V8.py 8 > analyze_V8.log 2>&1
../.venv/Scripts/python followup_V8.py > followup_V8.log 2>&1
FILES="certs_V8.jsonl.gz exp1_V8_certG.jsonl.gz"
[ -f exp2_V8_certG.jsonl.gz ] && FILES="$FILES exp2_V8_certG.jsonl.gz"
( time ../.venv/Scripts/python verify_joint.py 8 $FILES ) > verify_V8.log 2>&1
echo "exit=$?" >> verify_V8.log
echo CHAIN_DONE >> followup_V8.log
