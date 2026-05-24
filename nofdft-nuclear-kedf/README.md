# Nonlocal orbital-free kinetic-energy functionals for nuclei

Code, inputs and data for the study of a nonlocal, density-dependent-kF
kinetic-energy density functional (KEDF) for atomic nuclei, learned from
Kohn-Sham references with a quadratic tail-margin penalty and a selected-path
energy-rise penalty.

The point of this repository is that you can do the study yourself: build the
functional from the theory, run the calculations, get the same physics, and
then change something and find out what happens. It is not organized around
reproducing our output files byte for byte. Our numbers are recorded in
`PAPER_VALUES.md` as one run for orientation, and what counts as agreement is
stated in `tolerances.json`.

**Start here**: `RUNBOOK.md` walks the whole chain one calculation at a time.
`THEORY_TO_CODE.md` says what each piece of the model is and where it lives in
the code — read it first if you intend to reimplement rather than rerun.

## The physics in brief

The system is a simplified nuclear problem: N = Z, spherical, no spin-orbit,
no Coulomb, a simplified SkP interaction. The orbital-free kinetic-energy
density is written as a von Weizsacker term plus a nonlocal term built from a
dimensionless kernel evaluated at y = kF |r - r'| with a density-dependent
Fermi momentum, minus a local Thomas-Fermi-like correction. The kernel is the
Lindhard-based base kernel plus a learned correction delta_g(y); the local
correction carries one learned scalar.

The correction is fitted to Kohn-Sham references for O16, Ca40 and A=80 (A=140
is held out entirely) using pointwise kinetic data, the Euler-Lagrange
potential, an exact energy equality, and two penalties that are the substance
of the paper:

* a **quadratic tail-margin penalty**, which stops the correction from
  manufacturing spurious nonlocality-factor structure in the nuclear surface;
* a **selected-path energy-rise penalty**, which requires the energy to rise
  along explicitly constructed density directions. Without it, a shared
  one-kernel-for-all-nuclei fit looks good at the fit stage but is unstable:
  A=80 relaxes into a hollow, centrally depleted profile. That failure is a
  genuine fixed point — it survives a gradient regularization term and a finer
  angular quadrature — and only the fit-stage constraint removes it.

The main results: per-nucleus fits track the KS densities and nonlocality
factors closely; a single shared kernel does noticeably worse but transfers in
shape to the unseen A=140 system (structure and radius reproduced, absolute
energy overbound by a few percent).

## Layout

```
src/                     the solver (Fortran, one executable, namelist-driven)
namelists/               inputs by stage
  ks/                      Kohn-Sham reference runs
  teacher/                 base-kernel tables
  diagnostics/             value-only, diagonalization, Woods-Saxon starts
  preliminary/             the no-path correction diagnostics (templates)
production_inputs/       the eight production cases (ite / energy / nlf namelists)
inputs/ks/               KS densities used as OFDFT starting points
kernels/                 fitted kernels and their shifts
  per_nucleus/  universal/  diagnostics/
scripts/                 fits, tables, figures, helpers
  gp_data/                 base-kernel (teacher) tables and GP fit inputs
  path_data/               trial densities for the energy-rise penalty
results/                 our outputs
  raw/                     the eight production cases
  reference/               KS reference tables, KS energies, paper values
  preliminary_correction/  the no-path diagnostics
  tables/  figures/  gp/
tests/                   small consistency checks of the solver
reproduce.sh             optional wrapper for the production cases
```

Steps 1-4 of the runbook overwrite the shipped files under `inputs/`,
`scripts/gp_data/`, `scripts/path_data/` and `kernels/` with your own. That is
intended: it is what makes the chain use your results rather than ours. Keep
the downloaded archive, or work on a branch, if you want ours back.

## Requirements

gfortran with LAPACK/BLAS; python3 with numpy, scipy and matplotlib. No other
dependencies. A full pass through the chain is roughly 30-40 minutes of
single-core computation, dominated by the eight production cases.

## Running one thing, quickly

Build, then relax O16 with the released per-nucleus kernel:

```
cd src && ./build.sh && cd ..
mkdir -p work/quick && cd work/quick
cp ../../inputs/ks/ks_O16_nr250.dat           initial_density.dat
cp ../../kernels/per_nucleus/eg_PN14_O16.dat  kernel.dat
cp ../../production_inputs/per_nucleus/O16/ite.nml input_ks_nlf.dat
../../src/ks.x
```

`ite_density.dat` is the endpoint (column 2 is the species density; double it
for the total). It should sit essentially on top of the KS reference in
`../../inputs/ks/ks_O16_nr250.dat`. `RUNBOOK.md` STEP 5 explains the two
further readout stages that produce the published energy and NLF numbers.

## What is fixed and what is yours to change

Fixed by the model: the functional form, the base kernel, the exact energy
equality, and the two penalties. Everything else is a modelling choice, listed
with its value and its rationale in `THEORY_TO_CODE.md` (section 5) — basis
size and width, block weights, regularizers, the tail-margin target, the
energy-rise threshold, and the optional gradient coefficient. None of them
were tuned on the held-out A=140 system.

`THEORY_TO_CODE.md` section 7 gives the concrete recipes for adding a nucleus,
changing the Skyrme parametrization, replacing the basis or the penalties, and
for driving the released kernels from an independent orbital-free code.

## Honest limitations

* **Absolute energies do not transfer.** The shared functional overbinds the
  A=140 holdout by a few percent. Shape transfers; energy does not, at this
  basis size.
* **The A=80 shared solution is the most delicate result here.** Its central
  density moves by several percent between angular quadratures (nmu=24 to 48)
  and lands a few percent apart depending on whether the relaxation starts
  from the KS density or from a Woods-Saxon profile. Both are the same
  qualitative solution — centrally filled, no depletion — but do not read the
  third digit as physics.
* **The teacher tables shipped here came from an earlier run** of the same
  namelists and regenerate to about 1e-3 relative rather than exactly. That
  is the ordinary level of agreement for this quadrature, and it is the reason
  this repository judges agreement by tolerance rather than by file identity.
* **The fits are regularized least-squares problems with active sets.** Kernel
  coefficients and shift values are not observables; two visibly different
  kernels can give the same physics. Compare self-consistent results, not
  coefficients.
* **The energy-rise soft target is not met**, only the sign of the path energy
  is secured. That is sufficient for stability here but is not a proof of
  convexity along those directions.

## Development history in the repository

`results/preliminary_correction/` and `kernels/diagnostics/eg_NP14_shared.dat`
belong to the *preliminary* correction — the shared fit without the
energy-rise rows. They are kept because the instability they exhibit is the
reason the penalty exists, and because that claim should be checkable rather
than asserted. `namelists/preliminary/` holds the templates to regenerate
them, and RUNBOOK STEP 6 gives the commands.

## Citation and licence

See `CITATION.cff` and `LICENSE`.
