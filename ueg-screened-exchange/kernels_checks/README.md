# Executable checks (earlier kernels manuscript, partly used by the Mellin-pole paper)

These checks were written for an earlier manuscript on spectral kernels and Feynman-diagram weights, which was not published and whose material was revised into the Mellin-pole paper. The checks used by that paper are listed in `../mpp_checks/README.md` (SI Sec. S19). Section numbers below refer to the earlier manuscript.

All six scripts were executed for the final release on 29 September 2026.
Every stored numerical value agrees with the earlier archived record;
the normalization and independent-integral records now also contain an
explicit `PASS` status. The normalization, independent-integral and direct
Wick-enumeration scripts now enforce numerical or exact-count acceptance
checks rather than merely printing comparisons. Supplementary Material
Section S4 describes the results and their scope. The mapping of results to files is given below and in the repository README.

These files are self-contained and do not import the separate project code.
From the code archive root, run:

```sh
python -m pip install -r kernels_checks/requirements.txt
bash kernels_checks/run_all.sh
```

The supplied JSON files and text counts record the corresponding runs.
Running the scripts also produces detailed logs.
`environment.json` records the Python and package versions actually used:
Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0 and SymPy 1.14.0.
`requirements.txt` pins these four package versions. The scripts
can also be run individually and write results beside their source files.

| Script | Scope |
| --- | --- |
| `check_response_asymptotics.py` | Independent one- and two-dimensional self-energy quadratures; finite-part constant; large-width limit; local large-momentum frequency kernel; direct coupling/momentum integral for the leading low-density coefficient and its self-energy cancellation. |
| `check_exchange_wick_counts.py` | Direct Wick enumeration at interaction orders 2–4, including the two-chain topology at order 4. |
| `verify_exchange_chain_coefficients.py` | Independent enumeration of exchange and adjacent-chord topologies, with exact rational symmetry coefficients at orders 2–4. |
| `check_series_transfer.py` | Coupling-average and density-transfer identities on finite rational scalar kernels, together with power and logarithmic coupling moments. This is an algebra check, not an electron-gas energy calculation. |
| `normalization_checks.py` | Free-polarization normalization, the exchange-response energy prefactor, and order-by-order marking coefficients. |
| `independent_checks.py` | Frequency-kernel branches; leading low-density integral; self-energy finite parts; original occupied-to-empty transition coordinates versus radial reduction; the static small-momentum expansion; the high-density shell integral. |

The checks evaluate sampled identities and finite integrals. The uniform
spectral and density remainder estimates rely on the analytic arguments in
the manuscript. These scripts do not calculate the complete finite-density
RPA energy or establish the first surviving power of the combined response.
The starting diagram and its propagators remain those specified in the paper.

## Names in the archived verification code

The numerical integrals and their internal names have not been changed by
the final check. In the scripts, the Fock shapes `fx_float` and `fx_mp`
correspond to the manuscript's `f_L`; the scalar model quantities `e_kite_quad`
and `e_kite_closed` in the density-transfer check correspond to the coupling average denoted by `e_kite` in
the manuscript. The old term `recoil` in function names denotes the
large-momentum free-electron excitation scale. These are naming correspondences, not changes
to the integrals or normalizations tested by the scripts.
