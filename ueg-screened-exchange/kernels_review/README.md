# Reviewer verification notebook of the earlier kernels manuscript

This notebook was written for an earlier manuscript (Spectral kernels and Feynman-diagram weights for screened-exchange energies), which was not published and whose material was revised into the Mellin-pole paper. Appendix labels A–I refer to that manuscript. The notebook is kept as a record.

Execute the `.ipynb` cells in order.
The notebook is self-contained apart from its Python dependencies.

```sh
python -m pip install -r requirements.txt
python execute_notebook.py kernels_reviewer_verification.ipynb kernels_reviewer_rerun.ipynb
```

The executed notebook's final summary gives its assertion, failure and skip
counts. Tests cover
specified identities and sampled integrals, not formal proof of every analytic
estimate.

The screened line is `W_lambda`, including its numerator coupling, and
`Wbar = integral_0^1 W_lambda dlambda`; the energy is `2 K[Wbar]`.
The frequency variable is `nu`, the Lindhard shape is `L`, and the response
normalization is `Pi_0 = chi_0 = -q_tf_sq L`.
Vertex and self-energy response corrections share this unit, so their
logarithmic coefficients are `-q_tf_sq**2 alpha(z)/4` and
`+q_tf_sq**2 alpha(z)/4`. The geometric kernels and energy coefficients are
unchanged. The notebook distinguishes physical responses from their auxiliary
shapes and keeps source equation labels stable across the notation change.

The manuscript and bibliography
are not included here. The optional `mathematica/` route is a selected,
unexecuted Wolfram Language companion.
