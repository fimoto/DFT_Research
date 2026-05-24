# Mathematica checks for the JPA manuscript

This is an additional, Python-free review route for the revised JPA
manuscript. It checks selected principal formulas; it is not a full port of
the Python verification suite. No Wolfram kernel is available in the authoring environment:
the supplied files have **no claimed Mathematica execution results**.

Open `kernels_mathematica_review.nb` in Mathematica, start a fresh kernel, and
evaluate the input cell. All executable code is embedded in the notebook;
it needs neither Jupyter nor the companion `.wl` file. Alternatively evaluate:

```wolfram
Get["/absolute/path/kernels_checks.wl"]
```

Or from a terminal with Wolfram Engine configured:

```sh
wolframscript -file kernels_checks.wl
```

Inspect ``JPAReview`results`` and ``JPAReview`summary`` after execution.
Each check has a 120-second limit. A timeout or unresolved integral is reported
as TIMEOUT or UNRESOLVED, not as a successful comparison. Increase
`secondsPerCheck` in the source and rerun if necessary. NIntegrate messages
cause an unresolved result; numerical checks use scaled absolute error for
targets smaller than one. They are not rigorous interval enclosures.

## Correspondence and scope

The identifiers below are stable TeX source labels or named sections, so they
remain usable when the displayed equation numbers or pagination change.

| Group | Manuscript source labels / section | Calculation | Remaining analytic input |
|---|---|---|---|
| N | `app-eq:A0`, `app-eq:Onsager`; Fourier and Schwinger reduction | Normalization, Fourier transforms, regulated logarithm, D0 quadrature | Contraction; full coordinate reduction and integral interchange |
| J | `app-app-eq:Jclosed`, `app-eq:Wdef`, `app-eq:largeq-block` | Direct rational frequency integrals, both sign sectors; local profile `p_loc` | Occupation geometry |
| M | `app-app-eq:Ks`, `app-eq:largeMu`, `app-app-sec:layers` | Mellin kernel, Laurent coefficients, real moments | Weber-Schafheitlin expression is an input; profile endpoints and contour shift |
| L | `app-eq:A1-dispersion`, `app-eq:A1-beta`, `app-eq:A1-closed` | Dispersion integral substitutions, exact Beta integral, independent 2D quadrature | Full RPA to dispersion limit and uniform remainder |
| H | `app-eq:Az`, `app-eq:calA`, `app-eq:window`, `app-eq:Clog2` | Angular and frequency integrals, window and coupling integrals | Dominant-window proof |
| F | Exchange self-energy cancellation; `app-eq:fock-shell`, `app-eq:Edef` | Fock endpoints, shell integral, resolvent derivative, common-unit energy prefactor | Diagram signs and selected energy completion |
| S | `app-app-eq:staticavg`, `app-app-eq:c0`, `app-app-eq:c2`; static screening and density expansions | Generic Taylor integrand, exact coefficient integrals, parity, quantified bounds | Fixed-domain geometry, all-order analytic integration |
| B | `app-app-eq:lindhard-shell`, `app-app-eq:weighted-shell` | Exact endpoint mass and frequency-kernel bounds | Geometric integrated-shell remainder |

CAS results are conditional on explicit assumptions and on the equations chosen
as inputs. No comparison is evidence for those input equations themselves.
The full spectral reconstruction, finite-q Monte Carlo loop checks, spin
coefficient and the complete Python Check sections are not ported here. Consult the
existing equation coverage ledger for the broader Jupyter suite. The two
suites are complementary; this file is independently executable.

For publication, describe this as an **unexecuted Wolfram Language companion**
until a Mathematica user has run it and reviewed every non-PASS result.
This addition does not change any equation in the manuscript.

The source and its embedded notebook use the 2026-09-28 notation revision.
`qTFsq` denotes the squared reference Thomas--Fermi wavevector, and the local
integration variable `m` in J4 denotes the width ratio `s_Lambda`; its profile is
`p_loc`. F4 uses the common response normalization and the compensating energy
prefactor. The response rescaling does not change the kite energy.
Result labels identify current source labels or the mathematical quantity,
rather than carrying equation numbers from an earlier manuscript.

`VALIDATION.json` records the current source and notebook hashes, manuscript
fingerprints, and the exact-embedding check. This is a packaging check only;
it does not establish Wolfram syntax validity or evaluated CAS results.

The reporting harness uses [FullSimplify](https://reference.wolfram.com/language/ref/FullSimplify.html)
under explicit assumptions and [TimeConstrained](https://reference.wolfram.com/language/ref/TimeConstrained.html)
to retain inconclusive calculations without reporting success.
