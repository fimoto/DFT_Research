#!/bin/bash
# Regression test for the spatial-gradient regularization term.
# Checks its analytic EL potential against central finite differences with
# s2_grad active. Passes if relerr(TOTAL) < 1e-6.
set -e
BIN="${BIN:-../src/ks.x}"
[ -x "$BIN" ] || { echo "Build first: (cd ../src && ./build.sh)"; exit 1; }
cp ../inputs/ks/ks_O16_nr250.dat . 2>/dev/null || cp ../results/raw/ks_O16_nr250.dat . 2>/dev/null || true
cat > input_ks_nlf.dat << 'NML'
&input
  run_ks=.false. run_nlf=.false. run_derivative=.false. run_teacher=.false.
  run_deriv_check=.true. run_ite=.false. run_minimize=.false. run_diagscf=.false.
  nlf_use_ks_density=.false. external_density_file='ks_O16_nr250.dat'
  A=16 npart_q=8 kernel_mode='const' kf_const=0.80d0 s2_grad_in=0.5d0
  qmax=60.0d0 nq=6000 nmu=24 nr_tab=4001 use_semilocal=.false.
  deriv_check_file='dc_s2.dat'
/
NML
$BIN > dc_s2_log.txt 2>&1
err=$(grep -Eo "relerr\(TOTAL\)=[[:space:]]*[0-9.]+E[+-][0-9]+" dc_s2_log.txt \
      | grep -Eo "[0-9.]+E[+-][0-9]+" | sort -g | head -1)
echo "best relative gradient error (s2_grad=0.5): ${err:-see dc_s2_log.txt}"
awk -v e="${err:-1}" 'BEGIN{ if (e+0 < 1e-6){print "PASS";exit 0}else{print "CHECK";exit 1} }'
