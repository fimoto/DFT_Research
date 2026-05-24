#!/bin/bash
# Regression test: analytic kernel EL potential vs central finite differences.
# Passes if the derivative checker reports relative error below 1e-6.
set -e
BIN="${BIN:-../src/ks.x}"
[ -x "$BIN" ] || { echo "Build first: (cd ../src && ./build.sh)"; exit 1; }
cp ../inputs/ks/ks_O16_nr250.dat . 2>/dev/null || cp ../results/raw/ks_O16_nr250.dat . 2>/dev/null || true
cat > input_ks_nlf.dat << 'NML'
&input
  run_ks=.false. run_nlf=.false. run_derivative=.false.
  run_teacher=.false. run_deriv_check=.true. run_ite=.false.
  run_minimize=.false. run_diagscf=.false.
  nlf_use_ks_density=.false.
  external_density_file='ks_O16_nr250.dat'
  A=16 npart_q=8
  kernel_mode='const' kf_const=0.80d0
  qmax=60.0d0 nq=6000 nmu=24 nr_tab=4001
  use_semilocal=.false.
  deriv_check_file='deriv_check.dat'
/
NML
$BIN > deriv_check_log.txt 2>&1
# extract best (smallest) relerr(TOTAL) reported by the checker
err=$(grep -Eo "relerr\\(TOTAL\\)=[[:space:]]*[0-9.]+E[+-][0-9]+" deriv_check_log.txt \
      | grep -Eo "[0-9.]+E[+-][0-9]+" | sort -g | head -1)
echo "best relative gradient error: ${err:-see deriv_check_log.txt}"
awk -v e="${err:-1}" 'BEGIN{ if (e+0 < 1e-6) {print "PASS"; exit 0} else {print "CHECK deriv_check_log.txt"; exit 1} }'
