# Production solver inputs

This directory contains the exact namelist settings used by the repository's
production-calculation driver. They rerun the post-fit, self-consistent
orbital-free calculations from a released Kohn--Sham (KS) species density and
a fixed fitted correction kernel.  They do **not** refit the kernel or
regenerate the KS reference calculation.

`cases.tsv` is the compact index of all eight calculations. It gives the
source density, fitted kernel, physical case constants, and all three namelist
paths in one machine-readable row per case. The directories themselves are
the source of truth for the full solver configuration:

```text
production_inputs/
  per_nucleus/{O16,Ca40,A80,A140}/{ite,energy,nlf}.nml
  shared/{O16,Ca40,A80,A140}/{ite,energy,nlf}.nml
```

To inspect a fully staged case without compiling or running the solver, use:

```bash
./reproduce.sh --case per-nucleus/O16 --prepare-only
```

For a complete local execution and data-collection checklist, see
`RUNBOOK.md` in the repository root.

There are two kernel families and four symmetric nuclei:

| Family | Nucleus | `A` | `npart_q` | fitted correction file copied to `kernel.dat` | `c_corr_shift` | `s2_grad_in` | ITE monitor interval |
|---|---:|---:|---:|---|---:|---:|---:|
| `per_nucleus` | O16 | 16 | 8 | `kernels/per_nucleus/eg_PN14_O16.dat` | 0.020462 | 0.0 | 20 |
| `per_nucleus` | Ca40 | 40 | 20 | `kernels/per_nucleus/eg_PN14_Ca40.dat` | 0.083811 | 0.0 | 20 |
| `per_nucleus` | A80 | 80 | 40 | `kernels/per_nucleus/eg_PN14_A80.dat` | 0.080039 | 0.0 | 20 |
| `per_nucleus` | A140 | 140 | 70 | `kernels/per_nucleus/eg_PN14_A140.dat` | -1.198787 | 0.0 | 20 |
| `shared` | O16 | 16 | 8 | `kernels/universal/eg_SH14.dat` | -2.004040 | 0.0 | 50 |
| `shared` | Ca40 | 40 | 20 | `kernels/universal/eg_SH14.dat` | -2.004040 | 0.0 | 50 |
| `shared` | A80 | 80 | 40 | `kernels/universal/eg_SH14.dat` | -2.004040 | 0.0 | 50 |
| `shared` | A140 | 140 | 70 | `kernels/universal/eg_SH14.dat` | -2.004040 | 0.0 | 50 |

For each case, run the three stages in this order:

1. `ite.nml`: self-consistent imaginary-time evolution (ITE), initialized from
   `initial_density.dat` and writing `ite_density.dat`.
2. `energy.nml`: the higher-angular-quadrature (`nmu=64`) energy and functional-
   derivative evaluation on `ite_density.dat`.
3. `nlf.nml`: the `nmu=48` nucleon-localization-function evaluation on
   `ite_density.dat`.

Before running a case, copy
`inputs/ks/ks_<TAG>_nr250.dat` to `initial_density.dat` and the fitted file in
the table above to `kernel.dat` in an otherwise empty working directory.  The
executable always reads a file named `input_ks_nlf.dat`, so copy the desired
stage namelist to that name immediately before invoking `src/ks.x`.  Run every
stage in the same working directory.

The standardized primary outputs are:

| Stage | Primary outputs |
|---|---|
| ITE | `ite_history.dat`, `ite_density.dat`, `ite_summary.dat` |
| Energy/derivative | `energy.dat`, `derivative_radial.dat` |
| NLF | `nlf.dat` |

The solver also writes the diagnostic radial derivative and base-kernel tables
`ite_derivative_radial.dat`, `ite_kernel_radial.dat`, and
`derivative_kernel.dat`.  These are retained as useful numerical provenance;
they are not fitted-kernel inputs.

All cases use `kernel_mode='dd'`, `kf_const=0.80 fm^-1`, `qmax=60 fm^-1`,
`nq=6000`, `nr_tab=4001`, and no separately fitted semilocal term.  ITE uses
`nmu=24`, `ITE_2A0`, 8000 attempted steps, initial step
`1.0e-5 MeV^-1`, and adaptive-step parameters
`dt_max=5.0e-5`, `dt_grow=1.1`, and `dt_shrink=0.5`.  Setting `ite_tol=0`
deliberately disables early termination, so all 8000 attempted steps are made.

The production namelists also state the following values explicitly rather
than inheriting them from Fortran defaults:

| Input key | Value | Role |
|---|---:|---|
| `density_scale` | `1.0` | multiplicative scale applied while reading a density |
| `dd_xqmax_in` | `80.0` | dimensionless momentum cutoff used to build the density-dependent kernel table |
| `dd_ny_base_in` | `8001` | base number of points in that dimensionless table |
| `dd_kf_scale_in` | `1.0` | scale factor for the local density-dependent Fermi wave number |
| `kernel_lambda` | `1.0` | multiplier of the nonlocal kernel contribution |

Every file also sets `extra_kernel_file=''` and
`extra_g_file='kernel.dat'`: no separate dimensional-kernel correction is
loaded, and the staged dimensionless fitted correction is the only added
kernel shape.

The explicit `dd_*` keys control the density-dependent table used by ITE,
energy/derivative, and NLF stages. The separate legacy `teacher_*` keys are
used only when `run_teacher=.true.` to generate fitting-reference data; that
task is disabled in every production namelist here.
