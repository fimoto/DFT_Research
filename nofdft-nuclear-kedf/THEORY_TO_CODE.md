# From the theory to this code

This file is for someone who wants to understand what is being computed, or to
reimplement it. It states each ingredient of the model, then points at where
it lives in the code. Section numbers refer to the manuscript.

## 1. The energy functional

For N = Z nuclei without spin-orbit and without Coulomb, with
rho = rho_n + rho_p = 2 rho_q, the total energy is

    E[rho] = 2 * (hbar^2/2m) * Integral tau_q[rho_q] + E_Skyrme[rho]

with the simplified SkP interaction energy density

    0.375 t0 rho^2 + (t3/16) rho^(2+gamma) + (CG2/64) (grad rho)^2,
    CG2 = 9 t1 - 5 t2 - 4 t2 x2.

*Solver*: the interaction terms and their potentials are formed in the Skyrme
potential routine of `src/nofdft_ks_spherical_skp_nlf_v8.f90`; the constants
are set at the top of the file. *Python*: the same expressions appear in
`scripts/ks_total_energy.py` (energy) and in `kernel_fit.py` (`base_energy`,
`base_potential`), which is what lets the fit impose energy constraints
without calling the solver.

## 2. The kinetic-energy functional

The orbital-free kinetic-energy density is built as

    tau = tau_vW + tau_NL[g] - c_CORR * rho^(5/3)

where tau_vW is the von Weizsacker term, the nonlocal term is a double-radial
integral of a dimensionless kernel g evaluated at y = kF(r,r') |r - r'| with a
density-dependent Fermi momentum kF, and the last term is a local correction
whose coefficient is shifted by the fitted scalar s:

    c_CORR -> c_CORR + s.

The base kernel g_base comes from the Lindhard function; the fit adds
delta_g(y), tabulated in the released `kernels/**/*.dat` files as two columns
(y, delta_g).

*Solver*: selected by `kernel_mode='dd'`; the kernel table is read from
`extra_g_file` and the shift from `c_corr_shift`; `nmu` sets the angular
quadrature of the double integral, `dd_xqmax_in`/`dd_ny_base_in` the internal
tabulation of the base kernel. *Python*: `kernel_fit.py` reimplements the same
response — `base_dimensionless_kernel` (Lindhard), `tau_response` (the tau
produced by a given kernel on a given density) and `v_response` (its
functional derivative, including the rearrangement term from the density
dependence of kF).

The nonlocality factor plotted throughout the paper is

    NLF = 1 / (1 + chi^2),   chi = (tau - tau_vW) / tau_TF,

a bounded rewriting of the ratio the fit works with. Column 7 of the solver's
`nlf.dat` is the species NLF.

## 3. Solving for the density

Two independent solvers are provided, and agreement between them is what
establishes that a result is a fixed point rather than an artifact.

* `run_ite=.true.` — imaginary-time relaxation with an adaptive step size,
  rollback on energy increase, and particle-number restoration each step.
  Controlled by `ite_*` keys; writes density, history and summary files.
* `run_diagscf=.true.` — rearranged diagonalization: the Euler-Lagrange
  equation is written as a Schrodinger-like problem for sqrt(rho) with the
  nonlocal pieces treated as a potential, and iterated with mixing
  (`diag_lambda`, `diag_alpha`).

Other switches: `run_ks` (Kohn-Sham reference), `run_teacher` (base-kernel
tables), `run_nlf` (NLF readout on a fixed density), `run_derivative` (energy
and functional-derivative readout on a fixed density), `run_deriv_check`
(finite-difference test of the derivative).

The quantity monitored as convergence is the Euler-Lagrange residual R_EL, the
deviation of the functional derivative from a constant chemical potential
weighted by the density.

## 4. What the fit does

`scripts/kernel_fit.py` determines delta_g (14 coefficients on a basis of
TF-neutralized displaced Gaussians) together with the scalar shift s, by
solving a regularized least-squares problem with these blocks:

| block | what it says | code |
|---|---|---|
| value | reproduce chi from the KS reference on the KS density | `Xc`, `wch` in `load_nucleus` |
| EL potential | reproduce the functional derivative, up to one chemical potential per nucleus | `V`, `wEL`, `BETA` |
| exact energy | the total kinetic energy on the KS density is matched exactly | bordered KKT system in `solve` |
| surface band | weight the tail region separately | `mT`, `wT`, `ETAT` |
| tail margin | penalize chi dropping below 5 in the diagnostic band, so the correction does not create spurious NLF structure in the tail | active set in `solve`, `CHI_MIN` |
| energy rise | require the energy to rise along the trial paths of RUNBOOK STEP 3 (shared family only) | `path_rows`, `EPS_S` |

The two penalty blocks are the substantive modelling choices of the paper: the
tail-margin penalty is quadratic and one-sided, and the energy-rise rows are
what make the shared functional stable. Everything is solved by repeated
penalized least squares with a damped update (60 iterations, damping 1/2);
the active sets are re-evaluated each iteration.

The basis is TF-neutralized: each Gaussian has the Thomas-Fermi limit
projected out, so no member of the basis can shift the uniform-matter limit,
which is instead controlled solely by s.

## 5. Free parameters, and where they came from

| parameter | value | how chosen |
|---|---|---|
| basis size B | 14 | smallest that fits the KS chi profiles without oscillation |
| Gaussian width / centres | 0.9, y = 0..13 | resolve the kernel over the range of y reached in these nuclei |
| damping, iterations | 1/2, 60 | stability of the active-set iteration |
| ridge, shift regularizer | per-nucleus 3e-3, 3e-2; shared 1e-2, 1e-2 | condition the fit; exact energy equalities are enforced through the KKT border and are unaffected |
| EL weight eta | 2 (light), 12 (heavy) | balance the EL block against the value block |
| tail band weight | 1.2 (O16), 0.3 (others) | keep the surface region from being swamped |
| chi_min | 5 | soft target of the tail-margin penalty |
| energy-rise epsilon, weight | 2 MeV, 10 | the target is not reached and need not be; the sign of the path energy is what matters |
| gradient coefficient s2 | 0 in production | implemented and released, but it neither repairs the unstable preliminary fit nor improves the released one |

None of these were tuned on the A=140 system, which is held out from fitting
and from model selection.

## 6. Data formats

All data files are plain text with a comment header.

* KS species density (`inputs/ks/ks_<TAG>_nr250.dat`): two columns,
  r [fm] and rho_species [fm^-3], on r_i = i * rbox/(nr+1), rbox = 20 fm.
  Multiply by 2 for the total density.
* KS reference table (`results/reference/ks_nlf_<TAG>.dat`): r, rho_species,
  tau_TF, tau_vW, tau_KS, NLF_KS on the nr=500 grid.
* Teacher table (`scripts/gp_data/teacher_ddkf_<TAG>.dat`): nine columns; the
  fit uses column 1 (r) and column 7 (base density-dependent-kF tau).
* Kernel (`kernels/**/*.dat`): y and delta_g(y) on y in [0, 25]; the
  companion `_shift.txt` holds s as a single number.
* Solver NLF output (`nlf.dat`): column 1 is r, column 7 the species NLF.
* Endpoint density (`ite_density.dat`): r and species density, same
  convention as the input density.

If you supply your own Kohn-Sham reference, the first two formats are the only
interface you need to satisfy.

## 7. Extending the study

* **Another nucleus**: add a KS namelist under `namelists/ks/` with the right
  A and `npart_q`, run RUNBOOK STEPS 1-2, then add the tag to `NUCS` in
  `kernel_fit.py` (and to `BANDS`, `BETA`, `ETAT`, whose entries define the
  surface band and block weights for that mass).
* **Another Skyrme parametrization**: the constants live at the top of the
  solver source and are mirrored in `scripts/ks_total_energy.py` and
  `kernel_fit.py`. All three must agree, or the energy equality block will be
  fitting to a different functional than the solver integrates.
* **A different basis or penalty**: `tf_neutralized_basis` and the block
  assembly in `solve` are self-contained; the rest of the pipeline only sees
  the resulting kernel table and shift.
* **Your own solver**: the released kernels plus shift are a complete
  specification of the functional. Reading `kernels/universal/eg_SH14.dat`
  and its shift into an independent orbital-free code, and reproducing the
  qualitative statements at the end of `RUNBOOK.md`, is the strongest
  possible check on this work.
