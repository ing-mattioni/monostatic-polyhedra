#!/bin/sh
# calibration (15 -> 14) and spindle seeds at 13 / 14 ; all 3D objective G/vol^(1/3) with aspect cap 40
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=../.venv/Scripts/python
$PY search.py improve 15 --restarts 6 --maxevals 800000 --sigma 0.003 --norm 0 > logs_improve15.txt 2>&1
$PY search.py ladder --seedfile results/improve15_best.npy --restarts 2 --sigmas 0.003,0.03 --maxevals 400000 --norm 0 --amax 40 --tag lad15to14 > logs_lad15to14.txt 2>&1
$PY spindle.py 13 --n 30 --maxevals 400000 --tag spindle13 > logs_spindle13.txt 2>&1
$PY spindle.py 14 --n 18 --maxevals 400000 --tag spindle14 > logs_spindle14.txt 2>&1
$PY search.py cap --seedfile results/best14_G0.001939.npy --restarts 24 --maxevals 400000 --tag cap13 > logs_cap13.txt 2>&1
