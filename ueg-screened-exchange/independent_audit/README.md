# Independent re-derivation

`independent_audit.ipynb` is a second implementation of the quantitative content of both
manuscripts, written without importing `src/ueg_exchange`, `production/`,
`spin/reviewer_numerics.py` or any other module of this package, and without copying values
from the two reviewer notebooks. Every quantity is rebuilt from the definitions printed in
the manuscripts and only then compared.

It was produced by an AI assistant acting as an adversarial referee at the author's
request, and executed in full before distribution; `independent_audit.html` holds that
execution. The distributed run reports all items in agreement; the final cell raises if any disagree.

The notebook re-evaluates the kite energy at three densities. The ten-density
comparison quoted in the Supplemental Material came from the scripts in `grid/`,
which use the same estimator at larger sample sizes; see `grid/README.md` for
their settings and recorded output, and for what they do and do not retain.

It is supplementary evidence, not part of the locked production chain. It uses ordinary
floating point and modest sample sizes; runtime depends on the machine and BLAS threading.
Its Monte-Carlo errors are larger than those of `production/`. Sections 1-6
and 8-9 are deterministic; Section 7 is stochastic with fixed seeds.

What it covers, and what the reviewer notebooks cannot supply on their own:

- the factor chain of `A_0`, `D_0` by two routes, `e_2x`, the block endpoint coefficients,
  the Mellin kernel and the large-Lambda coefficients from the real-integral moments;
- `I_0`, `I_2`, `int C dq`, `C_inf` from the overlap-volume representation, hence `A_2`
  exactly;
- `F`, `J_C`, and a coupling-first reduction of `Phi_4` and `Psi` derived here and checked
  against the direct double quadrature, hence all three pieces of `A_3`;
- `S_3` shown to equal the (1856/675, 264/125) split identically in zeta, and the
  full-polarization rescaling applied to the three-term expansion;
- the static 4D reduction identity `A^2 - 4RS = 4(x+y)^2 H` as an exact polynomial
  identity, and the sum rule `int X_c dq = D_0`;
- an independent Monte-Carlo estimator whose bare limit must return the
  Onsager-Mittag-Stephen constant, evaluated at metallic densities and at `r_s = 1e8`
  against the three-term expansion;
- the `O(Q^-5)` cross-sector remainder of spin-paper SM Sec. S17 (formerly Appendix C), tested with the subtraction
  inside the integrand;
- the transcription of the ten published comparison values.

Requirements: `numpy`, `scipy`, `sympy`, `mpmath`. No pinned runtime; unlike the spin-paper
production entry points this notebook does not gate on interpreter version.
