# RUNBOOK — running the study yourself

The aim of this repository is not that you reproduce our output files
byte for byte. It is that you can carry out the same study independently:
build the same functional from the same theory, run the same calculations,
obtain the same physics, and then change things and see what happens.

So this runbook is written as a chain you actually walk. Each step reads what
the previous step produced — not a frozen copy of what we produced — and each
step tells you what a *successful* result looks like, as a range rather than a
digit string. The numbers we obtained are recorded separately in
`PAPER_VALUES.md`; they are one run on one machine, useful for orientation and
not a target to hit exactly. What agreement means quantitatively is declared
once, in `tolerances.json`.

The repository ships a complete set of starting-point files (KS densities,
teacher tables, trial-path densities, fitted kernels, results, tables,
figures) so that you can enter the chain at any step. Steps 1–4 overwrite
those files with your own; that is intended. Keep the downloaded archive (or
a `git stash`/branch) if you want to look back at ours.

Requirements: gfortran + LAPACK/BLAS for the solver; python3 with numpy,
scipy, matplotlib for the fits, tables and figures. Times below are for a
single ~2020s CPU core.

--------------------------------------------------------------------------
## STEP 0 — build the solver

```
cd src && ./build.sh && cd ..
```

Produces `src/ks.x`. All solver stages below are the same executable driven by
a namelist copied to `input_ks_nlf.dat` in the working directory; which task
runs is selected by the `run_*` switches. `THEORY_TO_CODE.md` maps those
switches to the equations.

--------------------------------------------------------------------------
## STEP 1 — Kohn–Sham reference solutions

The KS solutions play two roles: their densities are the starting points and
the fixed densities of the fits, and their kinetic-energy densities define the
target that the orbital-free functional is asked to reproduce.

For each TAG in O16 Ca40 A80 A140 (~50 s each):

```
mkdir -p work/ks_O16 && cd work/ks_O16
cp ../../namelists/ks/ks_O16_nr500.nml input_ks_nlf.dat
../../src/ks.x
python3 ../../scripts/resample_ks_density.py ks_O16_species.dat ks_O16_nr250.dat
cd ../..
```

Outputs: `ks_O16_species.dat` (nr=500 species density), `ks_nlf_O16.dat`
(the tau/NLF reference table), `ks_O16_levels.dat`, `ks_O16_history.dat`, and
the resampled `ks_O16_nr250.dat` used as the OFDFT starting density.

**Feed them into the chain** (this is what makes the rest of the runbook use
*your* results):

```
cp work/ks_O16/ks_O16_nr250.dat  inputs/ks/
cp work/ks_O16/ks_nlf_O16.dat    results/reference/
```

Then the consistent KS total energies:

```
python3 scripts/ks_total_energy.py
```

Expect roughly −126.5, −398.9, −895.7, −1677.4 MeV for O16, Ca40, A80, A140.
These are evaluated with the solver's own mass constant
ħ²/2m = 20.73553 MeV fm², so if you use a different constant they will move by
about 0.2%.

*If you implement the KS reference in your own code*, this is the natural
place to substitute it: everything downstream needs only the two files above,
whose formats are documented in `README.md` (Data formats).

--------------------------------------------------------------------------
## STEP 2 — base-kernel (teacher) tables

The fit does not learn the kinetic-energy density from scratch; it learns a
correction to a base nonlocal kernel. This step evaluates the base kernel on
the KS density and tabulates it, so the fit can work with differences.

For each TAG (~40 s each), in the directory holding `ks_<TAG>_species.dat`:

```
cd work/ks_O16
cp ../../namelists/teacher/teacher_O16.nml input_ks_nlf.dat
../../src/ks.x
cp teacher_ddkf_O16.dat ../../scripts/gp_data/
cd ../..
```

Column 7 is the density-dependent-kF base kinetic-energy density; the other
columns carry the constant-kF variant and the TF/vW pieces used by the
diagnostics.

A note on our shipped copies: they came from an earlier run of the same
namelists, and regenerating them here reproduces them to about 1e-3 relative
in column 7, not exactly. That is the ordinary level of agreement between
runs of this quadrature, and it is why the checks in this runbook are stated
as tolerances.

--------------------------------------------------------------------------
## STEP 3 — trial densities for the energy-rise constraint

**What this is for.** The shared (one-kernel-for-all-nuclei) fit needs one
extra ingredient. A shared fit built only from pointwise data reproduces the
KS kinetic-energy density well at the KS densities, yet the resulting
functional can still be *unstable*: relaxing A=80 from the KS density runs
downhill into a hollow, centrally depleted profile. The failure is invisible
at the fit stage because nothing in the pointwise data says anything about
energies away from the KS density.

The remedy is to give the fit a few explicit statements of the form "moving in
this direction must cost energy". This step constructs those directions:

```
python3 scripts/make_paths.py
```

It reads `inputs/ks/ks_A80_nr250.dat` and the recorded depleted endpoint
`scripts/path_data/a80_artifact_endpoint.dat`, and writes three trial
densities

    rho_xi = rho_KS + xi * (rho_depleted - rho_KS),   xi = 0.35, 0.60, 1.00

each nonnegative and carrying the same particle number as the KS density
(printed as `N_species` ≈ 40 for each). In STEP 4 the shared fit adds one row
per trial density requiring the total energy there to lie *above* the energy
at the KS density.

If you would rather derive the unstable direction yourself instead of using
our recorded endpoint, run the no-path fit first (STEP 4, family
`shared_nopath`), relax A=80 with it (STEP 6, preliminary diagnostics), and
point `make_paths.py` at your own endpoint density: any density that shows the
instability defines a usable direction.

Per-nucleus fits do not use these files at all.

--------------------------------------------------------------------------
## STEP 4 — fitting the kernel correction

```
python3 scripts/kernel_fit.py                    # every family, ~3 min
python3 scripts/kernel_fit.py pernucleus         # or one family at a time
python3 scripts/kernel_fit.py shared
python3 scripts/kernel_fit.py shared_nopath      # the no-path comparison fit
python3 scripts/kernel_fit.py no_tail value_only # diagnostic variants
```

**What each family needs.** Every family assembles the response of all four
nuclei, so all twelve files below are required even when you fit one nucleus;
only the shared family additionally uses the trial densities.

| family | needs | writes |
|---|---|---|
| `pernucleus` | `inputs/ks/ks_<TAG>_nr250.dat`, `results/reference/ks_nlf_<TAG>.dat`, `scripts/gp_data/teacher_ddkf_<TAG>.dat` (×4 each) | `kernels/per_nucleus/eg_PN14_<TAG>.dat` + `_shift.txt` |
| `shared` | the same twelve **plus** `scripts/path_data/path_A80_xi*.dat` | `kernels/universal/eg_SH14.dat` + `_shift.txt` |
| `shared_nopath` | the same twelve (no trial densities) | `kernels/diagnostics/eg_NP14_shared.dat` |
| `no_tail`, `value_only` | the same twelve | `kernels/diagnostics/eg_{NT14_*,VO14_O16}.dat` |

Each kernel is written together with a `_shift.txt` holding the scalar
local-correction shift; the solver needs both (`extra_g_file` and
`c_corr_shift`).

**What a good fit looks like.** `scripts/fit_report.txt` (rewritten on each
run, so use no arguments if you want all families in one file) reports, per
nucleus, the fixed-density NLF distance, the smallest chi in the surface
diagnostic band, and — for the shared family — the energy change along each
trial path. A fit is behaving if

* the per-nucleus fixed-density distances are a few times 1e-2 or below,
* the surface-band chi stays at or above the soft target of 5, so the
  correction does not manufacture spurious NLF structure in the tail,
* the shared trial-path energy changes are **positive** (ours are a few tenths
  of an MeV; the soft target of 2 MeV is not reached, and need not be — what
  matters is the sign).

The shift values and distances we obtained are in `PAPER_VALUES.md`. Yours
will differ in the last digits, and by more if your STEP 1–2 results differ.
The fit is a regularized least-squares problem with an active set; small input
changes move the coefficients smoothly, so treat these as continuous
quantities, not as fingerprints.

--------------------------------------------------------------------------
## STEP 5 — production self-consistent calculations

Eight cases: two families (per-nucleus, shared) × four nuclei. Each case is
three solver invocations in one directory. Running them by hand is worthwhile
at least once — you see exactly which density, which kernel and which angular
quadrature enters each published number.

**(a) Stage the case.** The namelists refer to the fixed names `kernel.dat`
and `initial_density.dat`, so the staging is two copies:

```
mkdir -p work/per_nucleus_O16 && cd work/per_nucleus_O16
cp ../../inputs/ks/ks_O16_nr250.dat           initial_density.dat
cp ../../kernels/per_nucleus/eg_PN14_O16.dat  kernel.dat
cp ../../production_inputs/per_nucleus/O16/ite.nml    .
cp ../../production_inputs/per_nucleus/O16/energy.nml .
cp ../../production_inputs/per_nucleus/O16/nlf.nml    .
```

For the shared family the kernel is `kernels/universal/eg_SH14.dat` and the
namelists come from `production_inputs/shared/<TAG>/`.

**Important if you refitted in STEP 4:** the shift lives in the namelists as
`c_corr_shift`, so after a refit you must copy your new value from
`kernels/.../<name>_shift.txt` into the three namelists of that case (or edit
`production_inputs/` once). `scripts/validate_production_inputs.py` checks the
namelists against `production_inputs/cases.tsv` and will flag inconsistencies
within the case set, but it cannot know which kernel you meant.

**(b) Stage 1 — the self-consistent calculation** (imaginary-time relaxation,
8000 steps, adaptive step size, nmu=24; ~90 s):

```
cp ite.nml input_ks_nlf.dat && ../../src/ks.x
```

Outputs: `ite_density.dat` (the endpoint density; column 2 is the species
density, so the central total density is 2 × column 2 at the first radial
point), `ite_history.dat` (energy, Euler–Lagrange residual, step size and
rollback flag at each monitor point), `ite_summary.dat`.

**(c) Stage 2 — energy readout** at the finer quadrature nmu=64 on the fixed
endpoint density (~15 s):

```
cp energy.nml input_ks_nlf.dat && ../../src/ks.x
```

Output `energy.dat`; its `E_total_NeqZ_MeV` line is the reported energy. The
stage exists because the energy monitored during stage 1 at nmu=24 is not
converged in the angular quadrature.

**(d) Stage 3 — NLF readout** at nmu=48 on the same fixed density (~10 s):

```
cp nlf.nml input_ks_nlf.dat && ../../src/ks.x
```

Output `nlf.dat`; column 7 is the species nonlocality factor.

**(e) Collect** under the names the table and figure scripts expect:

```
per-nucleus: ite_density -> results/raw/dens_pernuc_<TAG>.dat
             ite_history -> results/raw/hist_pernuc_<TAG>.dat
             energy      -> results/raw/pnE_<TAG>.dat
             nlf         -> results/raw/nlf_pernuc_<TAG>.dat
shared:      ite_density -> results/raw/v5_<TAG>_d.dat
             ite_history -> results/raw/v5_<TAG>_h.dat
             energy      -> results/raw/v5e64_<TAG>.dat
             nlf         -> results/raw/v5nlf_<TAG>.dat
```

`./reproduce.sh --case per-nucleus/O16 --output runs` performs exactly this
sequence if you prefer not to type it eight times; `--list` shows the cases
and `--prepare-only` stages without running.

**What a good result looks like.** Per nucleus, in the per-nucleus family, the
endpoint should sit on top of the KS reference: central density within a
percent or so, rms radius within ~0.01 fm, NLF distance of order 1e-2, and the
final Euler–Lagrange residual well below 0.1 MeV. In the shared family the
agreement is deliberately looser (one kernel for four nuclei): central
densities within roughly 10%, radii within a few 1e-2 fm, NLF distance around
1e-1 for the training nuclei. The qualitative statements that must survive any
correct implementation are listed at the end of this file.

--------------------------------------------------------------------------
## STEP 6 — diagnostics

These are the calculations behind Table I, the figures and the discussion,
rather than behind Tables II and III.

**Solver benchmark: the rearranged diagonalization** (Table I; ~2 min total).
This one depends on no fit at all — it uses the base kernel at constant kF —
so you can run it as soon as STEP 1 has given you a KS density. It is what
establishes that the Euler-Lagrange equation can be solved by a second,
independent method, and how the rearrangement parameter lambda and the mixing
alpha affect convergence.

```
mkdir -p work/diagbench && cd work/diagbench
cp ../../inputs/ks/ks_O16_nr250.dat .
for CFG in "1.0 1.0 l1_a1" "1.0 0.3 l1_a03" "3.0 1.0 l3_a1" "10.0 1.0 l10_a1"; do
  set -- $CFG
  sed -e "s/diag_lambda=3.0d0/diag_lambda=${1}d0/" \
      -e "s/diag_alpha=1.0d0/diag_alpha=${2}d0/" \
      -e "s/diag_l3_a1/diag_${3}/g" \
      ../../namelists/input_diag_benchmark.dat > input_ks_nlf.dat
  ../../src/ks.x
done
cd ../..
```

The fixed-step imaginary-time runs that the table compares against use the
same namelist with `run_diagscf=.false.`, `run_ite=.true.`,
`ite_dt_adapt=.false.`, `ite_dt=1.0d-5`, `ite_max_iter=8000` and `20000`, and
the output names `ite_<N>_{h,d,s}.dat`.

Collect the results into `results/raw/diag_benchmark/` and read the table off
them:

```
python3 scripts/make_diag_benchmark.py --raw results/raw/diag_benchmark
```

What the benchmark should show, for any correct implementation: the naive map
(lambda, alpha) = (1, 1) does **not** converge within 500 iterations, while
mixing at lambda = 1 or rearrangement at lambda = 3 or 10 restores
convergence; all converged settings reach the same energy (about
-117.135 MeV here) from the same starting density; and the fixed-step
relaxation approaches the diagonalization endpoint as the step budget grows
(tens of keV and a relative density difference of order 1e-2 at 8000 steps,
falling by roughly two orders of magnitude at 20000). The iteration counts and
the auxiliary-Hamiltonian gaps themselves are implementation details — ours
are in `PAPER_VALUES.md` — and the paper's point is only that the gap does not
by itself explain which settings converge.


**Value-only endpoint** (~2 min) — shows that fitting the pointwise kinetic
data alone gives a functional whose self-consistent solution drifts, here to a
central density well above KS:

```
mkdir -p work/vo && cd work/vo
cp ../../inputs/ks/ks_O16_nr250.dat .
cp ../../kernels/diagnostics/eg_VO14_O16.dat kernel.dat
cp ../../namelists/diagnostics/value_only_O16_ite.nml input_ks_nlf.dat && ../../src/ks.x
cp ../../namelists/diagnostics/value_only_O16_nlf.nml input_ks_nlf.dat && ../../src/ks.x
cd ../..
```

**Diagonalization cross-check** (~10 s per nucleus) — confirms that the
imaginary-time endpoints are genuine fixed points of the Euler–Lagrange
equation and not artifacts of the relaxation algorithm. Stage as above with
`kernels/universal/eg_SH14.dat` and `namelists/diagnostics/diagU_<TAG>.nml`.
Expect convergence (residual below 1e-8) to the same density as the relaxation
endpoint, agreeing to a relative L2 difference of order 1e-4 or better, with
the same central density to the digits you would print.

**Initialization robustness** (~10 s) — same check started from a Woods–Saxon
profile unrelated to the KS solution, using
`namelists/diagnostics/diagU_{O16,A80}_ws_init.nml`. For the light system the
two starts land on the same solution; for A=80 they land on the same
*qualitative* solution (centrally filled) but a few percent apart in central
density. That sensitivity is real and is discussed in the paper.

**Fixed-density diagnostic curves** (seconds):

```
python3 scripts/make_diagnostic_curves.py
```

**The preliminary (no-path) correction** (~10 min) — the demonstration that
the energy-rise rows of STEP 3 are necessary. With
`kernels/diagnostics/eg_NP14_shared.dat` staged as `kernel.dat`:

```
for S2 in 0.0 0.2 0.5 1.0 2.0; do
  sed -e "s/@S2@/$S2/" -e "s/@NMU@/24/" \
    namelists/preliminary/prelim_A80_ite_TEMPLATE.nml > input_ks_nlf.dat
  ../../src/ks.x
done
sed -e "s/@S2@/0.0/" namelists/preliminary/prelim_A80_diag_TEMPLATE.nml \
  > input_ks_nlf.dat && ../../src/ks.x
```

The expected finding is qualitative and robust: A=80 relaxes to a strongly
depleted centre (central density roughly a third of the KS value), the
depletion survives every value of the gradient coefficient s2 and also a finer
angular quadrature (`@NMU@=48`), and the diagonalization converges to the same
depleted density — i.e. it is a real fixed point of that functional, not a
relaxation artifact. Only adding the trial-path rows removes it.

--------------------------------------------------------------------------
## STEP 7 — tables, figures, and comparison with the paper

```
python3 scripts/make_tables.py --raw results/raw --ref results/reference --out results/tables
python3 scripts/gp_fit.py
python3 scripts/make_release_figures.py --gp-dir scripts --output-dir results/figures
```

(`gp_fit.py` first: the Gaussian-process figure reads its products.)

To see how your run compares with ours:

```
python3 scripts/validate_tables.py --tables results/tables \
    --expected results/reference/paper_example_values.json
```

Read its output as a comparison, not as a pass/fail gate. It reports the
differences entry by entry; judge them against `tolerances.json`. Differences
at or below those levels mean you have the same physics. Differences far
beyond them mean something substantive differs — most often the kernel used in
STEP 5 not matching the kernel fitted in STEP 4, or a shift value not copied
into the namelists.

--------------------------------------------------------------------------
## What should agree, and what need not

Quantities that are properties of the *theory* and should come out of any
correct implementation, at the level declared in `tolerances.json`:

* per-nucleus fits track the KS density and NLF closely, while a single shared
  kernel does noticeably worse on all four nuclei;
* the shared functional transfers in shape but not in absolute energy: the
  A=140 system, which enters no fit and no model selection, is reproduced in
  NLF structure and radius while being overbound by a few percent;
* a shared fit without the energy-rise rows is unstable at A=80, the
  instability is a genuine fixed point, and neither a gradient term nor a
  finer quadrature repairs it;
* fitting pointwise kinetic data alone yields a functional whose
  self-consistent density drifts away from KS;
* dropping the surface-tail data and the tail-margin penalty lets spurious NLF
  structure appear in the surface region.

Quantities that are properties of *our particular run* and should not be
expected to match beyond the tolerances: individual kernel coefficients and
shift values, energies at the 0.1 MeV level, central densities in the third
decimal, and everything in `PAPER_VALUES.md`.
