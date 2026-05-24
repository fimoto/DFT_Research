# pySecDec check of two coefficient integrals

This directory uses an external, unrelated code base, pySecDec (sector
decomposition for parametric integrals, https://github.com/gudrunhe/secdec),
to reproduce two closed-form coefficient integrals of the spin paper with
an *automated* treatment of their endpoint behaviour. Nothing from `src/`,
`production/` or the reviewer notebooks enters.

| package    | integral | closed form | pySecDec (recorded run) |
|------------|----------|-------------|-------------------------|
| `fp_cross` | cross-sector finite part `Fp int_0^inf x^-4 J_C(x^-4) dx`, spin-paper Eq. (A3cross-profile) | `-(8/5) sqrt(2 pi) \mathcal L(3/4)/\mathcal L(1/4) = -1.35554093567` | `-1.35554093601 +- 9.6e-10`, `delta^-1` coefficient `0` |
| `me_euler` | Euler integral `M_E` of the Lindhard piece, spin-paper Appendix A | `(176/75) sqrt(2 pi) \mathcal L(3/4)/\mathcal L(1/4) = 1.98812670565` | `1.98812670560 +- 8.1e-11` |

Why these two and not the others. Sector decomposition acts on integrands that
are products of polynomials raised to powers on the unit hypercube, and it
resolves singularities on coordinate faces. The finite part above is exactly
that: after `x = t/(1-t)` and the coupling-integral form of `J_C`, the integrand
is `2 l t^(-2+delta) (1-t)^(2-delta) (t^4 + l (1-t)^4)^(-1/2)`, and the hand
subtraction of the `4/(3 eps)` power divergence in the manuscript is replaced by
the analytic regulator `delta`: pySecDec must return a vanishing `delta^-1`
coefficient and the finite part as the `delta^0` term. That is an independent,
mechanised confirmation of the subtraction. `M_E` is a Beta-type integral with
an integrable `w^-1/2` endpoint and needs no subtraction; it is included as a
second, cheap point of contact with the same tool.

The six-dimensional kite integral itself is not a sector-decomposition problem:
its domain is the symmetric difference of two balls, its Coulomb singularity
lies on the interior diagonal `p = k`, and its screened line contains the
non-polynomial Lindhard function. Rewriting it into hypercube polynomial form
would reproduce the analytic reduction by hand, after which pySecDec adds
nothing over its own integrators. The kite energy is therefore checked by the
three stochastic implementations described in the Supplemental Material, not
here.

## Reproduce

    pip install pySecDec            # tested with pySecDec 1.6.6, Python 3.12, g++ 13
    python3 generate.py             # writes fp_cross/ and me_euler/ (generated C++)
    sh make_all.sh                  # compiles both (a few minutes)
    python3 run.py                  # evaluates and compares; prints ALL OK

The generated package directories are large and are not distributed; only the
generator, the driver, and `recorded_output.txt` are. The Qmc integrator's
error is a spread over random lattice shifts and is approximate; the recorded
deviations are 1-2 of those errors.
