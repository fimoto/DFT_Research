# Numerical checks for the Mellin-pole paper

*Revisiting the $r_s$ expansion of the uniform electron gas in $d$ dimensions: $\log r_s$ terms from coincident Mellin poles* (submitted). Manuscript version v4.21: main text Eqs. (1)–(59), Supporting Information (SI) Eqs. (S1)–(S104).

These scripts reproduce the numbers of the paper that come from quadrature, series evaluation or fits, rather than from a closed form. The closed forms and the symbolic identities are checked by the Wolfram Language notebook in `../mpp_review/`. The checks of the three-dimensional kite that were written for the earlier kernels manuscript are in `../kernels_checks/` and are listed in the table below.

Each script writes `results/<script>.json`. That file lists every reproduced number with the place where it is printed, the printed value, the value computed here and a verdict. For a printed decimal the check passes if the printed digits are the correctly rounded or the truncated value. For a printed bound, such as "agrees to $10^{-17}$", it passes if the deviation found here does not exceed the bound.

```sh
python -m pip install -r mpp_checks/requirements.txt
python mpp_checks/run_all.py             # all scripts (about 1.5 h on one core)
python mpp_checks/run_all.py --quick     # skips c06, c07 and c12, small Monte Carlo (about 3 min)
python mpp_checks/run_all.py --full-mc   # Monte Carlo with 4 x 1e8 samples, as in SI Sec. S19
python mpp_checks/fig_data.py            # data behind main-text Figs. 1, 2 and 4
```

The stored results were produced with `--full-mc` under Python 3.13, NumPy 2.3.5, SciPy 1.17.0, mpmath 1.3.0 and SymPy 1.14.0. Each JSON file records the versions that were used. The scripts import only these packages and `common.py`. They do not depend on `src/`.

## Where each result is checked

| Manuscript item | Script | Notes |
| --- | --- | --- |
| SI Table S2; SI Sec. S1: $h_d$ from the hypergeometric form vs the angular average ($d=2,3,3.5,4,5$, 18 digits); main Eq. (10); SI Eqs. (S3), (S4); $g_2(0.3)$, $g_2(1)$, $g_2(3)$ and the vanishing $q$, $q^3$ coefficients (40 digits); $g_3$ of main Sec. 2.9 | `c01_lindhard_d.py` | |
| Main Eqs. (30)–(32), Table 2 (all columns); SI Sec. S2: Eq. (S5), $(d-3)A_d/2$ at $d=3.1,3.01,3.001$; $A_d/d$ at $d=100,1000,3000$; $-1/(2e\pi)$; Eq. (21) of Schlesier et al. for $3.001\le d\le3000$; closed form of $\Sigma_D$ to 40 digits | `c02_ring_d.py` | Uses Eqs. (18), (21) and (B10) of Schlesier et al., Phys. Rev. B 102, 035123 (2020) |
| Main Sec. 2.9, Eqs. (36)–(37); SI Sec. S3: $R_3$ quadrature at $z=0.1,0.3,1,3$ (printed values, $2.5\times10^{-17}$); $\int hR_3\,dz$, $\int h^2g_3\,dz$; $C$ to 60 digits; SI Sec. S19: the 50-digit quadratures of main Eqs. (30) and (37) | `c03_response_3d.py` | |
| Main Eq. (38) and properties (i)–(iii) of Sec. 3.1; SI Sec. S4: Fourier ingredients, large-$z$ series, angular quadrature at $z=0.3,1,3$ ($7\times10^{-18}$, $10^{-13}$), Ward slopes $0.318246$, $0.1429308$, $-0.01930873$; main Eq. (58) and the coefficient $0.0389$ of Table 4 | `c04_R2.py` | about 2 min |
| Main Eqs. (39), (40), (59); SI Sec. S5: $\mathsf A$, $\mathsf B$, $\mathsf C$, quadratures at 40, 60 and 100 digits, PSLQ, value in Ry and Hartree, size at $r_s=0.1$ | `c05_a2l.py` | |
| SI Sec. S6, Table S3: $R(q,z)-R_2(z)$ at $z=0.3$, $1$ and $q=0.08,0.04,0.02$; the factor about $3.6$ on halving $q$; the sign change at $z=1$ | `c06_finiteq_2d.py`, `finiteq2d.py` | about 3 min. Independent integrator (polar coordinates about $\bm p$, elementary radial integral). This integrator needs $\nu=qz>0$; the static values are in `c12`. |
| SI Sec. S7: $[R_2(1;\kappa_{\rm s})-R_2(1)]/\kappa_{\rm s}$ for $\kappa_{\rm s}=10^{-1}$–$10^{-5}$ and $-0.0321699$ | `c07_screened_line.py` | about 7 min |
| SI Sec. S8, Table S4: both rows, the fully polarized values, the small-$\zeta$ coefficients, the limit $0.81649$ and its approach | `c08_spin_2d.py` | |
| SI Sec. S19: Lindhard closed form vs three-dimensional quadrature at $(0.7,0.4)$ and $(3,2)$, 30 digits; $A_0$; Monte Carlo of $\mathcal K[v]$; main Eq. (43) | `c09_kite_3d.py` | The Monte Carlo uses its own seed. It reproduces the sampling scheme and the statistical agreement with $e_{2x}$, not the printed digits $0.04836(2)$ |
| Main Sec. 3.3, SI Sec. S19, SI Table S5: table entries from the stored unrounded values; free fits $c_{\rm fit}=-0.00091(5)$, $-0.00092(12)$; residual sums $0.02$ and $380$; $c_{1,\rm fit}$ with $c_{\rm fit}=c_{\rm hd}$ | `c10_hd_fit.py` | Reads `validation/high_density/results/` |
| SI Sec. S6, static comparison with Table I of Das Sarma: $\Delta R/R(0)$ at $x=0.1$–$0.4$, printed $4.523\times10^{-4}$, $2.104\times10^{-3}$, $5.141\times10^{-3}$, $9.699\times10^{-3}$ | `c12_static_dassarma.py`, `finiteq2d_static_tips.py`, `finiteq2d_static.py` | about 55 min. `finiteq2d_static_tips.py` ends the $\phi$ panels at the crescent tips $\cos\phi_0=-q/2$, where $E_{\bm p}=0$, grades the nodes towards them and extrapolates in the number of nodes. `finiteq2d_static.py` is the earlier code with the panels split at $\pm\pi/2$. It is not converged at $z=0$ and is kept as the record of the row printed in v4.20 |
| Main Sec. 3.3, Eq. (53), Table 3: $A_1$ from the $(\lambda,x)$ integral, $q_*^{-1}$, the Table 3 coefficients, $A_1/q_*$ at $r_s=10^8$, $q_*k_F$ at $r_s=2\times10^5$, $10^6$; the stored direct values at $r_s=10^8$, $2\times10^5$, $10^6$ and the BRM values quoted in the text | `c11_lowdensity_kite.py` | Reads `production/raw/tableI/` and `production/reference_brm.json` |
| SI Sec. S19: closed forms of $J(\omega;a,b)$ ($2.2\times10^{-16}$); the double radial integral (S75) by two methods ($1.5\times10^{-14}$) and in the original coordinates (orders 48–128); the large-momentum self-energy kernel at $w_\Lambda=0.5,2,10$, $q=8,32,128$; $p_{\rm loc}$; $A_1=-0.21469376098909768$, the moment $0.25000000000000006$ and $-A_1$; the coupling identities ($8.5\times10^{-12}$) | `../kernels_checks/` (`independent_checks.py`, `check_response_asymptotics.py`, `check_series_transfer.py`, `normalization_checks.py`) | Stored in `../kernels_checks/*.json` |
| Main Sec. 3.3 and SI Table S5: the direct values of $e_{\rm kite}$ at $r_s=10^{-3}$–$10^{-1}$; $c_1=0.01150$~Ry | `../validation/high_density/run_direct_high_density.py`, `hd_coefficients.py` | Companion paper (spin paper) |
| Main Sec. 3.3: the momentum cutoff near $12k_F$ | `../validation/run_cutoff_scan.py`, `analyze_cutoff_scan.py` | Companion paper |
| Main Figs. 1, 2 and 4: the plotted quantities | `fig_data.py` | Writes `figure_data/*.csv`. The figures are drawn from these files by `../mpp_figures/` |

## What these checks do not establish

The scripts evaluate the integrals and series that the paper reports as numerical evidence. They do not prove the remainder bounds, SI Propositions S1 and S2 or SI Lemma S1, and they do not test hypotheses (H1) and (H2). The finite-$q$ values of `c06` and the slopes of `c04` and `c07` are numerical observations at the stated points, as in the paper.
