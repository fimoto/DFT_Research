#!/bin/bash
# Regression test: detect omission of the spatial-gradient regularization
# energy from the total-energy summary. Evaluates the SAME density with
# s2_grad=0 and s2_grad=0.5, and checks that the change in E_total_NeqZ_MeV
# equals an independent discrete integral of the two-species gradient energy.
set -e
BIN="${BIN:-../src/ks.x}"
[ -x "$BIN" ] || { echo "Build first: (cd ../src && ./build.sh)"; exit 1; }
cp ../inputs/ks/ks_O16_nr250.dat . 2>/dev/null || cp ../results/raw/ks_O16_nr250.dat . 2>/dev/null || true
run() {  # s2 value -> writes egy_$1.dat
cat > input_ks_nlf.dat << NML
&input
  run_ks=.false. run_nlf=.false. run_derivative=.true. run_teacher=.false.
  run_deriv_check=.false. run_ite=.false. run_minimize=.false. run_diagscf=.false.
  nlf_use_ks_density=.false. external_density_file='ks_O16_nr250.dat'
  A=16 npart_q=8 kernel_mode='const' kf_const=0.80d0 s2_grad_in=$1
  qmax=60.0d0 nq=6000 nmu=24 nr_tab=4001 use_semilocal=.false.
  derivative_energy_file='egy_$2.dat'
/
NML
$BIN > /dev/null 2>&1
}
run 0.0d0 0
run 0.5d0 5
python3 - << 'PY'
import re, numpy as np
def get(f,key): return float(re.search(key+r'\s+(-?\d+\.\d+E[+-]\d+)',open(f).read()).group(1))
E0=get('egy_0.dat','E_total_NeqZ_MeV'); E5=get('egy_5.dat','E_total_NeqZ_MeV')
Eg_summary=get('egy_5.dat','E_grad_species_MeV')  # per species
# independent discrete integral of the two-species gradient energy on the SAME density
h2m=20.73553; s2=0.5
d=np.loadtxt('ks_O16_nr250.dat'); r=d[:,0]; dr=r[1]-r[0]; rho=d[:,1]
rp=np.zeros_like(rho); rp[1:-1]=(rho[2:]-rho[:-2])/(2*dr); rp[0]=(rho[1]-rho[0])/dr; rp[-1]=(rho[-1]-rho[-2])/dr
W=4*np.pi*r**2*dr
Eg_indep_2sp = 2.0*h2m*s2*np.sum(W*np.maximum(rho,1e-20)**(1/3)*rp**2)
dE = E5-E0
print("E_total(s2=0)   = %.6f"%E0)
print("E_total(s2=0.5) = %.6f"%E5)
print("Delta E_total   = %.6f"%dE)
print("2*E_grad(summary)= %.6f"%(2*Eg_summary))
print("independent 2sp integral = %.6f"%Eg_indep_2sp)
ok1 = abs(dE - Eg_indep_2sp) < 1e-4
ok2 = abs(2*Eg_summary - Eg_indep_2sp) < 1e-4
# energy closure: reconstruct E_total from printed components
Tvw=get('egy_5.dat','T_vW_species_MeV'); Tnl=get('egy_5.dat','T_KEDFnl_species_MeV')
Tsl=get('egy_5.dat','T_SL_fit_species_MeV'); Eskp=get('egy_5.dat','E_SkP_total_MeV')
Erec = 2*(Tnl + Tsl + Eg_summary) + Eskp
ok3 = abs(Erec - E5) < 1e-6
print("closure: 2(T_KEDFnl+T_SL+E_grad)+E_SkP = %.6f vs E_total %.6f"%(Erec,E5))
print("PASS" if (ok1 and ok2 and ok3) else "FAIL")
import sys; sys.exit(0 if (ok1 and ok2 and ok3) else 1)
PY
