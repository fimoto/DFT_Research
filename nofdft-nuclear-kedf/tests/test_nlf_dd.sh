#!/usr/bin/env bash
# Regression test for the density-dependent-kF NLF postprocessor.
# Verifies:
#   (1) run_nlf with kernel_mode='dd' + nlf_use_ks_density=.false. uses the dd
#       convolution (not the constant-kF fallback);
#   (2) the NLF read from run_nlf equals the NLF written by run_derivative on the
#       same density (both use the SCF kernel table and quadrature) to a tight tolerance;
#   (3) NLF is finite and in [0,1] everywhere.
# Requires: a built solver (../src -> ks.x) and a converged density file $DENS.
set -euo pipefail
BIN="${BIN:-./ks.x}"
DENS="${DENS:-dens_pernuc_O16.dat}"
EG="${EG:-eg_PN14_O16.dat}"
SHIFT="${SHIFT:-0.020462}"
A="${A:-16}"; NQ="${NQ:-8}"

# (1)+(3) run_nlf dd
cat > _t_nlf.dat << NML
&input
  run_ks=.false. run_nlf=.true. run_derivative=.false. run_ite=.false.
  nlf_use_ks_density=.false. external_density_file='$DENS'
  A=$A npart_q=$NQ kernel_mode='dd' kf_const=0.80d0
  extra_g_file='$EG' c_corr_shift=${SHIFT}d0 s2_grad_in=0.0d0
  qmax=60.0d0 nq=6000 nmu=48 nr_tab=4001 nlf_file='_t_nlf_out.dat'
/
NML
cp _t_nlf.dat input_ks_nlf.dat
OUT=$("$BIN" 2>&1)
echo "$OUT" | grep -q "dd NLF output" || { echo "FAIL: run_nlf did not use the dd path"; exit 1; }
echo "$OUT" | grep -q "constant-kF convolution matrix cache" && { echo "FAIL: dd run fell back to constant-kF"; exit 1; }

# (2) run_derivative on the same density (writes NLF in radial cols 9,10)
cat > _t_der.dat << NML
&input
  run_ks=.false. run_nlf=.false. run_derivative=.true. run_ite=.false.
  nlf_use_ks_density=.false. external_density_file='$DENS'
  A=$A npart_q=$NQ kernel_mode='dd' kf_const=0.80d0
  extra_g_file='$EG' c_corr_shift=${SHIFT}d0 s2_grad_in=0.0d0
  qmax=60.0d0 nq=6000 nmu=48 nr_tab=4001
  derivative_energy_file='_t_der_e.dat' derivative_file='_t_der_r.dat'
/
NML
cp _t_der.dat input_ks_nlf.dat
"$BIN" > /dev/null 2>&1

python3 - << 'PY'
import numpy as np, sys
nlf=np.loadtxt('_t_nlf_out.dat')          # col 7 (idx6) = NLF
der=np.loadtxt('_t_der_r.dat')            # col10 (idx9) = NLF
r=nlf[:,0]; a=nlf[:,6]; b=np.interp(r,der[:,0],der[:,9])
if not np.all(np.isfinite(a)):            sys.exit("FAIL: NLF not finite")
if a.min()<-1e-9 or a.max()>1+1e-9:       sys.exit("FAIL: NLF out of [0,1]: [%.3f,%.3f]"%(a.min(),a.max()))
d=np.max(np.abs(a-b))
print("max|NLF_run_nlf - NLF_derivative| = %.3e"%d)
if d>1e-6:                                sys.exit("FAIL: run_nlf and run_derivative NLF disagree (%.3e)"%d)
print("PASS: dd NLF path consistent with SCF derivative path; NLF in [0,1].")
PY
rm -f _t_nlf.dat _t_der.dat _t_nlf_out.dat _t_der_r.dat _t_der_e.dat input_ks_nlf.dat
