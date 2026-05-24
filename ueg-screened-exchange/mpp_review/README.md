# Wolfram Language checks for the Mellin-pole paper

`MPP_checks_v4.21.ipynb` checks the closed forms and the symbolic integral identities of *Revisiting the $r_s$ expansion of the uniform electron gas in $d$ dimensions: $\log r_s$ terms from coincident Mellin poles* (manuscript version v4.21: main text Eqs. (1)–(59), SI Eqs. (S1)–(S104)). `MPP_checks_v4.21.wls` holds the same code cells as a script.

Every check is symbolic (`Simplify`, `Integrate`, `Limit`, `Series`, `Residue`, `Solve`) or evaluates a closed form to the printed digits. The notebook does not use `NIntegrate`. Numbers that the paper obtains only by quadrature, series evaluation or fits are checked by the Python scripts in `../mpp_checks/`. Each check is named by the printed equation, table or section number of v4.21. The cells follow the order of the derivations, so a cell uses only results checked above it.

## Running

1. Install the free Wolfram Engine (free for personal and developer use) and activate it with `wolframscript`.
2. For Jupyter, install the Wolfram Language kernel (WolframLanguageForJupyter, `./configure-jupyter.wls add`), open the notebook, choose the "Wolfram Language" kernel and run all cells.
3. Without Jupyter, run `wolframscript -file MPP_checks_v4.21.wls` in this directory.

Both write `wolfram_check_results_v4.21.json` with one verdict per check. `exact`, `numerically equal (20 digits)`, `not evaluated (numerically equal)` and `matches printed` are passes. Every other verdict needs attention.

## Status

The notebook was synced from the v4.18 notebook. Every equation, table, figure, section and proposition number was mapped through the label maps of the two manuscript versions. The checks of content that v4.20 and v4.21 no longer contain were removed: the chain family and its energy expression (former Eq. (44), SI Eqs. (S45) and (S48), SI Table S6) and the display on the low-density weight of a high-frequency tail (former Eq. (57)). The remaining checks are unchanged, except that the check of the static comparison in SI Sec. S6 follows the converged values printed in v4.21.

**This version of the notebook has not been run yet.** The results file will be added after it has been executed.
