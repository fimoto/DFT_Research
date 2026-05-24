# The numbers in the paper — one run, for orientation

Everything below is what *our* run produced, on one machine with one compiler.
It is recorded so you can see whether your own run landed in the same place,
and so the figures and tables in the paper can be traced to specific files.
It is **not** a set of targets to reproduce digit by digit. What counts as the
same result is declared in `tolerances.json`; the qualitative statements that
should survive any correct implementation are listed at the end of
`RUNBOOK.md`.

Solver constants: hbar^2/2m = 20.73553 MeV fm^2; simplified SkP interaction
(t0 = -2931.70, t1 = 320.62, t2 = -337.41, t3 = 18708.97, x2 = -0.53732,
gamma = 1/6), N = Z, no spin-orbit, no Coulomb, spherical symmetry.

## Kohn-Sham references

| system | A | N_species | E_KS [MeV] | rho_c [fm^-3] | r_rms [fm] |
|---|---|---|---|---|---|
| O16  | 16  | 8  | -126.496  | 0.1448 | 2.741 |
| Ca40 | 40  | 20 | -398.893  | 0.1858 | 3.396 |
| A80  | 80  | 40 | -895.744  | 0.1769 | 4.092 |
| A140 | 140 | 70 | -1677.411 | 0.1976 | 4.807 |

## Fitted kernels

Local-correction shifts (`_shift.txt`, entering the namelists as
`c_corr_shift`), basis size B = 14 TF-neutralized displaced Gaussians:

| fit | shift | max abs delta_g |
|---|---|---|
| per-nucleus O16  | +0.020462 | 0.334 |
| per-nucleus Ca40 | +0.083811 | 0.119 |
| per-nucleus A80  | +0.080039 | 0.781 |
| per-nucleus A140 | -1.198787 | 0.919 |
| shared (with energy-rise rows) | -2.004040 | 1.213 |
| shared, no-path comparison fit  | -1.100868 | 0.806 |

Fit-report quantities: per-nucleus fixed-density NLF distances 0.0125 /
0.0206 / 0.0266 / 0.0292; shared 0.108 / 0.090 / 0.094 / 0.079 (O16 / Ca40 /
A80 / A140). Shared trial-path energy changes +0.41, +0.29, +1.09 MeV at
xi = 0.35, 0.60, 1.00 — all positive, which is the point.

## Table II — per-nucleus functionals

`E64` is evaluated at nmu=64 on the converged nmu=24 endpoint density;
parentheses give the KS reference.

| system | E_KS | E64 | abs(dE) | R_EL | rho_c (KS) | r_rms (KS) | d |
|---|---|---|---|---|---|---|---|
| O16  | -126.496  | -126.885  | 0.02 | 0.017 | 0.1459 (0.1448) | 2.741 (2.741) | 0.0090 |
| Ca40 | -398.893  | -400.543  | 0.34 | 0.054 | 0.1876 (0.1858) | 3.402 (3.396) | 0.0145 |
| A80  | -895.744  | -898.140  | 0.20 | 0.045 | 0.1812 (0.1769) | 4.104 (4.092) | 0.0248 |
| A140 | -1677.411 | -1683.405 | 4.08 | 0.049 | 0.2026 (0.1976) | 4.811 (4.807) | 0.0271 |

## Table III — shared functional (s2 = 0, no gradient term)

A140 enters no fit and no model selection.

| system | role | E64 | dE_KS | rho_c (KS) | r_rms (KS) | d |
|---|---|---|---|---|---|---|
| O16  | training | -127.844  | -1.3  | 0.1619 (0.1448) | 2.774 (2.741) | 0.0960 |
| Ca40 | training | -401.125  | -2.2  | 0.1719 (0.1858) | 3.404 (3.396) | 0.0965 |
| A80  | training | -898.913  | -3.2  | 0.1834 (0.1769) | 4.088 (4.092) | 0.0960 |
| A140 | holdout  | -1719.083 | -41.7 | 0.1903 (0.1976) | 4.793 (4.807) | 0.0739 |

Holdout overbinding: 41.7 MeV, 2.5% of the total.

## Table I — rearranged-diagonalization benchmark

Constant-kF O16, nr=250, kF=0.80 fm^-1, nmu=24; iteration counts for
R_rho < 1e-8 fm^-3/2; gap between the two lowest eigenvalues of the final
frozen auxiliary Hamiltonian.

| lambda | alpha_mix | iterations | E [MeV] | gap [MeV] |
|---|---|---|---|---|
| 1  | 1.0 | no convergence in 500 | — | — |
| 1  | 0.3 | 50  | -117.1351065 | 11.365 |
| 3  | 1.0 | 41  | -117.1351065 | 10.852 |
| 10 | 1.0 | 136 | -117.1351064 | 10.411 |

Fixed-step imaginary time against the (3, 1.0) diagonalization endpoint:
8000 steps give dE = 41.6 keV and a relative L2 density difference of
1.0e-2; 20000 steps give 0.37 keV and 4.6e-4.

## Diagnostics

* Diagonalization cross-check of the shared endpoints: converged with residual
  below 1e-8; relative L2 density difference to the relaxation endpoint
  2.2e-5 / 7.9e-5 / 1.9e-4 / 2.9e-4 (O16 / Ca40 / A80 / A140), same central
  density to four digits.
* Woods-Saxon initialization: O16 reaches the KS-initialized fixed point to
  1e-7; A80 reaches a centrally filled solution but a different one,
  rho_c = 0.197 versus 0.183 (relative L2 1.5e-2).
* Angular quadrature, shared A80: going from nmu=24 to nmu=48 moves rho_c from
  0.1834 to 0.1977 (7.8%), the radius by 0.004 fm, d from 0.096 to 0.094, and
  the nmu=64 energy readout by 0.12 MeV.
* Gradient term, shared A80: s2 = 0.5 fm^4 moves rho_c from 0.1834 to 0.1988,
  i.e. away from the KS value 0.1769. Hence s2 = 0 in production.
* Value-only fit, O16 endpoint: rho_c = 0.1860 against the KS 0.1448.
* No-tail fits, fixed-density NLF: spurious surface structure peaking at 0.60
  near 6.2 fm (O16) and 0.86 near 7.1 fm (Ca40).
* Preliminary (no-path) shared correction, A80 relaxed from KS: rho_c =
  0.0608 / 0.0609 / 0.0612 / 0.0617 / 0.0629 for s2 = 0, 0.2, 0.5, 1.0, 2.0
  fm^4; 0.0626 at nmu=48; the diagonalization converges to the same depleted
  fixed point (rho_c = 0.0608, relative L2 to the relaxation endpoint 2.0e-4).
* Gaussian-process reformulation: selected hyperparameters (1.5, 0.3, 300),
  shift -8.32492, maximum posterior-mean correction 0.098, 2-sigma coverage of
  the A140 holdout 0.12.

## Where these files live

`results/raw/` holds the eight production cases under the canonical names
listed in RUNBOOK STEP 5; `results/tables/` the generated tables;
`results/figures/` figures 1-6; `results/reference/` the KS reference tables,
the consistent KS energies, and `paper_example_values.json` (a machine-readable
copy of Tables II and III above); `results/preliminary_correction/` the no-path
diagnostics; `results/gp/` the Gaussian-process products.
