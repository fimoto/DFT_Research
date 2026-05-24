# ueg-screened-exchange

Code and data for two papers by F. Imoto. Directory names refer to the papers by content:

* **spin paper** (`spin/`, `plot_data/spin/`, `REPRODUCE_SPIN.md`): *Low- and high-density asymptotics of RPA-screened second-order exchange in the spin-polarized uniform electron gas* (submitted);
* **Mellin-pole paper** (`mpp_checks/`, `mpp_review/`, `mpp_figures/`, and the checks in `kernels_checks/` listed in its table below): *Revisiting the $r_s$ expansion of the uniform electron gas in $d$ dimensions: $\log r_s$ terms from coincident Mellin poles* (submitted), manuscript version v4.21.

`kernels_checks/`, `kernels_review/` and `plot_data/kernels/` were written for an earlier manuscript on spectral kernels and Feynman-diagram weights. That manuscript was not published. Its material was revised into the Mellin-pole paper, which uses some of these checks (SI Sec. S19). The rest is kept as a record and is listed in the section on the earlier manuscript below.

Every number, table and figure of these two papers and their supplementary material can be regenerated from this repository, either from stored raw data or by recomputation.

## Quick start

```bash
python3.12 -m venv .venv && . .venv/bin/activate
python -m pip install -r requirements-spin.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONHASHSEED=0
PYTHONPATH=src python tests/test_all.py     # unit and contract tests
python spin/analytic_audit.py                # exact checks of the analytic results
python validation/analyze_validation.py      # SM S8-S11 summary from the stored replicates
```

The production records were written with Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. `production/verify.py` and `production/compute.py` refuse other interpreter versions before writing results. Run scripts from the repository root (the `validation/` scripts change to it themselves).

## Spin paper: where each result comes from

| Manuscript item | Code | Data |
|---|---|---|
| Coefficients $A_1,A_2,A_3$, spin moments, Eqs. (31)–(58) | `spin/analytic_audit.py`, `spin/verify_spin_*.py`, `scripts/verify_analytics.py`, `scripts/check_A3_highprec.py` | `spin/generated/` |
| Tables I, III and IV, Figs. 4 and 5, SM Figs. S3–S5 and Table SXX (Estimator A, direct integration) | `production/compute.py`, `production/verify.py`, `production/render.py` | `production/raw/` |
| SM Sec. S15 and Table SXVI, fixed-grid quadrature | `production/compute.py` (fixed grids of SM Sec. S15; data files keep their historical `appendixB` names) | `data/appendixB_grid_*.dat` |
| Fig. 4 radial profile, SM S10 (Estimator B) | `spin/reviewer_numerics.py` | `figs/data/`, `validation/results/axis/` |
| Figs. 2 and 4–6, SM Figs. S2–S5 | `figs/make_*.py` | `figs/data/`, `plot_data/spin/` |
| SM S1–S6 (sampling, errors, external check) | `production/`, `external/`, `pysecdec_check/` (requires pySecDec) | `production/raw/`, `pysecdec_check/recorded_output.txt` |
| SM S7 (third implementation) | `independent_audit/grid/`, `validation/run_s7_converged.py` | `validation/results/s7_converged.json` |
| Table IV column B and SM Sec. S10 (Tables SVIII–SIX): Estimator B interpolation residuals on the converged grid | `validation/run_residual_grid.py`, `run_residual_grid_ext.py`, `analyze_residual_ext.py` | `validation/results/residual_grid_rs10.json`, `residual_grid_ext.json` |
| Sec. V C and SM Sec. S22 (Table SXXI): momentum-cutoff test of the difference from the published values | `validation/run_cutoff_scan.py`, `analyze_cutoff_scan.py` | `validation/results/cutoff_scan.json`, `cutoff_scan_summary.json` |
| SM S8–S11, seed check in Sec. V A | `validation/` (see `validation/README.md`) | `validation/results/` |
| Fig. 6 and Table II, high-density expansion, Secs. IV and V B, SM Sec. S19 and Tables SXVII–SXIX: $c_{\rm hd}(\zeta)$, $c_1(\zeta)$, direct evaluations at $r_s=10^{-3}$–$10^{-1}$ | `validation/high_density/hd_coefficients.py`, `run_direct_high_density.py`, `fit_high_density.py`, `make_fig_high_density.py` | `validation/high_density/results/` |

Figures 1 and 3 are diagrams drawn in the manuscript source. `REPRODUCE_SPIN.md` gives step-by-step commands for the production tables and figures.

## Earlier kernels manuscript (not published): material kept as a record

Figure and section numbers in this table refer to that manuscript and its Supplementary Material (SM).

| Item | Code | Data |
|---|---|---|
| Constants and coefficients ($e_{2x}$, $A_0$, $D_0$, $c_{\rm hd}$, $A_1$, Mellin and series coefficients) | `src/ueg_exchange/`, `scripts/verify_analytics.py`, `tests/test_all.py` | `data/coefficients.csv`, `data/block_kernel.csv` |
| Supplementary Material Sec. S4 (independent checks, normalization, response asymptotics, Wick counts, chain coefficients, series transfer) | `kernels_checks/run_all.sh` | `kernels_checks/*.json`, `kernels_checks/wick_counts.txt` |
| SM Fig. S1, profile $\Phi(y)$ | `figs/make_phi_fig.py` | `figs/phi_grid_corrected.npy` |
| SM Figs. S2 and S3, $S(\Lambda)$, pole dispersion and $p_{\rm loc}$ | `figs/make_lorentzian_figs.py` | `figs/phi_grid_corrected.npy`, `figs/data/` |
| SM Fig. S4, static kernel $X_c(q)$ | `figs/make_xc_fig.py` | `figs/xc_grid.npy` |
| Grid of $\Phi(y)$ | `scripts/sample_phi_grid.py`, `scripts/recompute_phi_range.py`, `scripts/recompute_phi_fig1_range.py`, `scripts/build_phi_grid_corrected.py`; checks `scripts/phi_point_crosscheck.py`, `scripts/check_phi_panels.py`, `scripts/phi_grid_point_checks.py` | `figs/phi_grid.npy`, `figs/phi_recomputed_*.json`, `figs/phi_fig1_range_recomputed.json`, `figs/phi_grid_point_checks.csv` |
| Grid of $X_c(q)$ | `scripts/sample_xc_grid.py` | `figs/xc_grid.npy` |
| End-to-end reduction and Appendix B checks | `scripts/check_reduction_endtoend.py`, `scripts/check_companion_appendixB.py`, `scripts/check_ring_no_ln2.py` | printed |
| Reviewer verification notebook, Mathematica cross-checks | `kernels_review/` | in the notebook |

Run times on one core: `bash kernels_checks/run_all.sh` about 10 s (it rewrites the stored records byte for byte, apart from `environment.json`, which records the local package versions); the reviewer notebook about 7 min (`python -m pip install -r kernels_review/requirements.txt`, then `python kernels_review/execute_notebook.py kernels_review/kernels_reviewer_verification.ipynb rerun.ipynb`, 388 checks); `scripts/build_phi_grid_corrected.py` a few seconds; the point checks of the $\Phi$ grid several minutes each; the Monte Carlo sampling of the grids hours.

## Mellin-pole paper: where each result comes from

Equation, table and section numbers refer to manuscript version v4.21 (main text Eqs. (1)–(59), Supporting Information (SI) Eqs. (S1)–(S104)). `mpp_checks/README.md` lists every reproduced number with its script, and each script stores its checks in `mpp_checks/results/`.

| Manuscript item | Code | Data |
|---|---|---|
| Closed forms and symbolic identities of the main text and SI | `mpp_review/MPP_checks_v4.21.ipynb` (Wolfram Language; `.wls` script version) | `mpp_review/wolfram_check_results_v4.21.json` after execution |
| $d$-dimensional constants and the Lindhard shape $h_d$, $g_2$, $g_3$ (main Sec. 2.2, SI Sec. S1, Table S2) | `mpp_checks/c01_lindhard_d.py` | `mpp_checks/results/` |
| Ring series in $d$ dimensions, $A_d$ and Table 2, comparison with Schlesier et al. (main Sec. 2.8, SI Sec. S2) | `mpp_checks/c02_ring_d.py` | `mpp_checks/results/` |
| First-order exchange polarization in 3D and the Carr–Maradudin coefficient (main Sec. 2.9, SI Sec. S3) | `mpp_checks/c03_response_3d.py` | `mpp_checks/results/` |
| Closed form of $R_2$, its quadrature check and the Ward slopes (main Sec. 3.1, Eq. (58), SI Sec. S4) | `mpp_checks/c04_R2.py` | `mpp_checks/results/` |
| Coefficient $a_{2\ell}$ (main Eqs. (39)–(40), SI Sec. S5) | `mpp_checks/c05_a2l.py` | `mpp_checks/results/` |
| Finite wave number in 2D, SI Table S3 and the static comparison with Das Sarma (SI Sec. S6) | `mpp_checks/c06_finiteq_2d.py`, `mpp_checks/finiteq2d.py`, `mpp_checks/c12_static_dassarma.py`, `mpp_checks/finiteq2d_static*.py` | `mpp_checks/results/` |
| Screened internal exchange line (SI Sec. S7); spin polarization, SI Table S4 (SI Sec. S8) | `mpp_checks/c07_screened_line.py`, `mpp_checks/c08_spin_2d.py` | `mpp_checks/results/` |
| Kite normalization, Lindhard function and Monte Carlo of $\mathcal K[v]$ (SI Sec. S19) | `mpp_checks/c09_kite_3d.py` | `mpp_checks/results/` |
| Frequency integrals, radial reductions, $A_1$ and its cancellation, coupling identities (SI Sec. S19) | `kernels_checks/` (see also the table of the earlier manuscript above) | `kernels_checks/*.json` |
| Direct high-density values of $e_{\rm kite}$ and the fits (main Sec. 3.3, SI Sec. S19, Table S5) | `validation/high_density/run_direct_high_density.py`, `mpp_checks/c10_hd_fit.py` | `validation/high_density/results/`, `mpp_checks/results/` |
| Low-density term $A_1/q_*$, Table 3, direct values at $r_s=10^8$, $2\times10^5$, $10^6$ and the momentum-cutoff comparison with BRM (main Sec. 3.3) | `mpp_checks/c11_lowdensity_kite.py`, `production/`, `validation/run_cutoff_scan.py` | `production/raw/tableI/`, `production/reference_brm.json`, `validation/results/cutoff_scan_summary.json` |
| Data behind main Figs. 1, 2 and 4 | `mpp_checks/fig_data.py` | `mpp_checks/figure_data/` |
| Main Figs. 1–4 (Figs. 1, 2, 4 matplotlib, Fig. 3 TikZ) | `mpp_figures/make_figures.py` | `mpp_figures/fig*.pdf` |

Run times on one core: `python mpp_checks/run_all.py` about 1.5 h (`--quick` about 3 min); the Wolfram notebook a few minutes.

The data behind every plotted curve of the spin paper and of the earlier kernels manuscript are in `plot_data/spin/` and `plot_data/kernels/`. `LIMITATIONS.md` states what each numerical check does and does not establish, including the sampling of the $\Phi$ grid.

## Integrity

`release_manifest.json` holds the SHA-256 hash of every file in the GitHub release (the files admitted by `.gitignore`):

```bash
python production/verify.py --manifest
```

## Layout

```
src/ueg_exchange/    core library (Lindhard function, kernels, constants)
production/          Estimator A (direct integration), verification, raw records
external/            integrand of the representation used by production/
spin/                spin-paper analytic checks and Estimator B (imaginary-axis quadrature)
independent_audit/   third, separately written implementation (SM S7)
validation/          spin-dependent and convergence checks (SM S7-S11)
scripts/             analytic checks, grid generation for Phi and X_c, diagnostics
figs/, plot_data/    figure scripts and the data behind every plotted curve
kernels_checks/      checks written for the earlier kernels manuscript, partly used by the Mellin-pole paper
kernels_review/      reviewer notebook and Mathematica checks of the earlier kernels manuscript
mpp_checks/          Mellin-pole paper: numerical checks, figure data
mpp_review/          Mellin-pole paper: Wolfram Language notebook of the closed forms
mpp_figures/         Mellin-pole paper: scripts of the main-text figures and the PDF files
pysecdec_check/      external check of two coefficient integrals
tests/               unit and contract tests
data/                input and reference tables
```

