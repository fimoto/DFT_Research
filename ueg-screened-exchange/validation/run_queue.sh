#!/bin/bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
cd "$(dirname "$0")/.."
V=validation/results
S=20261001
# 1. full-polarization identity: zeta=1 at rs, zeta=0 at rs*2^(-4/3), same scale and seed
python3 - <<'PY'
import math
for rs in [0.1,1,10,100,1e4]: print(rs, rs*2**(-4/3))
PY
for rs in 0.1 1 10 100 10000; do
  rs2=$(python3 -c "print(repr($rs*2**(-4/3)))")
  python3 validation/run_spin_validation.py --rs $rs --zeta 1 --m 19 --reps 16 --seed $S --scale 2 --out $V/identity
  python3 validation/run_spin_validation.py --rs $rs2 --zeta 0 --m 19 --reps 16 --seed $S --scale 2 --out $V/identity
done
# 2. Fig.8 convergence: nested m=17,19,21
for m in 17 19 21; do
  python3 validation/run_spin_validation.py --rs 1 10 60 --zeta 0 0.4 0.8 1 --m $m --reps 16 --seed $S --scale 2 --out $V/fig8_m$m
done
python3 validation/run_spin_validation.py --rs 10 --zeta 0 0.4 0.8 1 --m 19 --reps 32 --seed 20262001 --scale 2 --out $V/fig8_m19_seedB
# 5. low-density spin dependence
python3 validation/run_spin_validation.py --rs 1000 5000 20000 100000 --zeta 0 0.4 0.8 1 --m 19 --reps 16 --seed $S --scale qstar --out $V/lowdensity
# 2b. near-full polarization
python3 validation/run_spin_validation.py --rs 10 60 --zeta 0.95 0.99 --m 19 --reps 16 --seed $S --scale 2 --out $V/nearfull
echo QUEUE_DONE
