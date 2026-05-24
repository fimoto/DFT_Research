# Spin-paper reviewer notebook

From this directory run `python run_inprocess.py`, or restart the Jupyter
kernel and Run All. Use Python 3.12.14 and `../requirements-spin.txt`.
See `../REPRODUCE_SPIN.md` for separate numerical computation and rendering steps.

This code/data distribution intentionally excludes manuscript sources. The
notebook regenerates production TeX exports in temporary directories and reports
the manuscript-source hash check as skipped. Included generated TeX tables and
numerical macros are numeric inputs rather than manuscript source. Mathematical and numerical tests run
normally. `equation_map.json` identifies the aligned manuscript and its equation numbers.
`ANALYTIC_COVERAGE.md` maps the current numbered displays to the checks that
test them (definition, symbolic identity, numerical check or partial coverage). It is built
from the executed notebook by `python build_coverage_ledger.py`; `python verify_analytic_coverage.py`
rebuilds it and fails if it is stale, if a check fails, or if a printed label cites an equation
number that the manuscript does not have. Printed labels must start with `Eq.(n)` so that the
ledger can read them.

Check numbering follows the argument, not the cell order: Check 31 collects the
exact derivation certificates for Checks 1-27 and is therefore placed before the
slow independent integrators of Checks 28-30. The execution summary counts every
check regardless of position.


## Added spin-interpolation checks

The matching article and its worked derivation notes are supplied in the
separate manuscript archive. Section V includes the positive leading endpoint-interpolation
bias, its analytic sign proof and the limitation of density-independent
interpolation. From the code root, run:

```bash
python figs/make_spin_interpolation.py
python spin/verify_spin_interpolation.py
python figs/make_spin_bias.py
python spin/verify_spin_bias.py
python spin/verify_spin_bias_algebra.py
```

Figure 7 illustrates analytic coefficients. Figure 8 reanalyses 104 saved
spin-map records with paired, covariance-aware replicate standard errors;
it performs no new integral. Figure 9 shows scaled full-energy residuals.
The spin-bias CSV and generated numerical macros record the finite-density
examples. The standalone checks distinguish algebraic identities and the
continuous sign argument from stored-data and supplementary numerical checks.

## Common notation

The full screened interaction is `W_lambda = lambda*v/epsilon_lambda`;
its coupling average is `Wbar = v R(c)/2`, and `e_kite = 2 K[Wbar]` in
Ry per electron. `R(c)` is a dimensionless screening factor, not the
coupling-averaged interaction itself. The positive Lindhard shape `L` is
distinct from the density response `Pi_0 = -q_tf_sq L` and from the
one-electron Green function `G`. Imaginary frequency is `nu`.
The notebook's numerical estimator and historical production data retain
their original normalization and sample values; the physical response
definitions and their displayed notation match the revised article.
