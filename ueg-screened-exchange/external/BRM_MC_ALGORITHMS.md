# BRM direct Monte Carlo algorithms

This note documents the Monte Carlo estimators implemented for the
Benites--Rosado--Manousakis (BRM) kite formulas, Eqs. (28)--(32) of
Phys. Rev. B 110, 195151 (2024).  The goal is reproducibility and auditability:
each estimator targets the same coupling-integrated screened second-order
exchange correction,

```text
epsilon_2b = epsilon_2b^0 + Delta epsilon_2b,1 + Delta epsilon_2b,2,
```

where `epsilon_2b^0` is the Onsager bare exchange value.  Energies are reported
in Ry per electron.

The BRM paper specifies the post-frequency-integration formulas and the
compactification `nu = a tan u`, but it does not specify the random sequence,
the infinite-`q1` map, any finite cutoff, the proposal densities, or the batch
definition used for the quoted error bars.  The algorithms below therefore
distinguish the integral identity from the sampling diagnostic.

## Common conventions

After rotating `q1` onto the positive `z` axis, the nominal 11-dimensional BRM
integral becomes a 9-dimensional integral.  The factor `4*pi` from the
orientation of `q1` is included in the radial measure.

The standard variables are

| Coordinate | Meaning |
| --- | --- |
| `lambda` | coupling constant in `[0,1]` |
| `tq` | radial map variable for `q1 in [0,infinity)` or `[0,qmax]` |
| `k` | point in the unit Fermi ball, subject to `|k+q1|>1` |
| `q2` | point in the unit Fermi ball, subject to `|q2+q1|>1` |
| `u` | frequency variable, `nu=a tan(pi u/2)` |

For `q=q1` along `z`,

```text
a(q,p) = q^2/2 + q p_z.
```

The screened factor is BRM Eq. (32),

```text
f_lambda^sigma(q,nu,zeta)
 = x_sigma^3 S_sigma(q,nu)
   / [x_sigma^2 + lambda alpha rs S_sigma(q,nu)/(pi q^2)],
```

with

```text
S_sigma = sum_sigma' x_sigma' g((x_sigma/x_sigma') q,
                                (x_sigma/x_sigma')^2 nu).
```

The unpolarized and polarized cases use the same code path; a spin component
with `x_sigma=0` is skipped.

## Algorithm A: literal uniform-ball sampler

This is the closest reconstruction of the formula as printed.

1. Sample `k` uniformly in the unit ball.
2. Sample `q2` uniformly in the unit ball.
3. Reject the sample unless both BRM particle-hole conditions hold:

   ```text
   |k+q1| > 1,     |q2+q1| > 1.
   ```

4. Keep the Coulomb denominators explicitly:

   ```text
   Delta1 denominator: |k-q2|^2
   Delta2 denominator: |k+q1+q2|^2
   ```

This is useful as a diagnostic because it closely follows the published
formula.  It is not a good production estimator.  Sampling a three-dimensional
Coulomb singularity with a flat proposal produces heavy tails; ordinary
Gaussian standard errors can be misleading at accessible sample counts.

Python implementation:

```bash
python3 mc/brm_literal_mc.py --proposal literal --q-map rational
```

Simple Fortran implementation:

```bash
cd fortran
make simple_brm_kite_mc
./simple_brm_kite_mc 5.0 0.0 4096 8 2.0 1 2 20260807
```

The seventh argument `2` selects the literal sampler.

## Algorithm B: exit-cap and Coulomb-centered importance sampler

This estimator evaluates the same integral with much lower variance.

1. Sample `k` uniformly in the exact exit cap

   ```text
   C(q) = { k : |k| <= 1, |k+q zhat| > 1 }.
   ```

   Its volume is

   ```text
   Vcap(q) = 4*pi/3,                      q >= 2,
           = pi q (12 - q^2)/12,          0 < q < 2.
   ```

2. For `Delta1`, sample `q2` around the singular point `q2=k`:

   ```text
   q2 = k + rho n,      rho in [0,2].
   ```

3. For `Delta2`, sample `q2` around the singular point `q2=-k-q1`:

   ```text
   q2 = -k - q1 + rho n,   rho in [0, |-k-q1|+1].
   ```

4. Use a radial-uniform density in `rho`.  The factor
   `rho^2 d rho d Omega` cancels the explicit Coulomb denominator
   `1/rho^2`.

This is the preferred independent estimator for the full infinite-domain BRM
integral.

Python implementation:

```bash
python3 mc/benites_direct.py --rs 5 --zeta 0 --m 16 --replicates 8
python3 mc/brm_literal_mc.py --proposal importance --q-map rational
```

Simple Fortran implementation:

```bash
cd fortran
make simple_brm_kite_mc
./simple_brm_kite_mc 5.0 0.0 4096 8 2.0 1 1 20260807
```

The seventh argument `1` selects the importance sampler.

## Algorithm C: maps for the transfer momentum

The BRM paper does not state how the unbounded `q1` radial integral was mapped.
The comparison code therefore implements three maps:

| `q_map` | Formula | Meaning |
| --- | --- | --- |
| `1` / rational | `q = s t/(1-t)` | infinite integral |
| `2` / tangent | `q = s tan(pi t/2)` | infinite integral |
| `3` / finite | `q = qmax t` | diagnostic hard cutoff |

For the infinite maps, the exact integral is independent of the scale `s`.
Residual scale dependence at finite sample size is a convergence diagnostic.

The finite map is not an alternative representation of the same integral.  It
computes the truncated cumulative contribution up to `qmax`.  It is useful for
asking how much high-`q1` tail is needed to reproduce a published number.

Example:

```bash
./simple_brm_kite_mc 1000000.0 0.0 65536 16 12.0 3 1 20260807
```

This evaluates the importance-sampled estimator with a hard cutoff
`q1 <= 12 k_F`.

## Algorithm D: independent scrambles or pseudo-random replicates

The Python code uses randomized Sobol sequences. Each replicate uses an
independent scramble, and the two BRM terms use independent scrambles.
Reported uncertainties are empirical standard errors,
`sqrt((s1^2+s2^2)/R)`, where `st^2` is the sample variance of the `R`
replicate means for term `t`. For differences between map nodes sharing
the same scrambles, form the difference within each term and replicate
before computing its variance. Finiteness of the production sampler's
second moment has not been proved. These errors quantify between-replicate variability,
not rigorous confidence intervals or bounds on systematic error.
The simple Fortran code uses the compiler's `random_number`
pseudo-random generator, with independent seeds for each term and replicate.

The Fortran code is intentionally not a bitwise reproduction of the Python
Sobol estimator.  It is a transparent scalar reference implementation.  For
high-precision production checks, prefer the Python randomized-Sobol estimator
or the existing MPI Fortran code.

## Recommended convergence checks

For any quoted number, check at least:

1. sample count per replicate;
2. number of independent replicates;
3. `q_scale` for infinite maps;
4. rational versus tangent map;
5. finite-cutoff cumulative values, especially at large `r_s`;
6. literal versus importance sampler on small test cases;
7. spin endpoints `zeta=0` and `zeta=1`.

At large `r_s`, the relevant moving shell is

```text
q_* = (4 q_TF^2 / 3)^(1/4).
```

The finite-cutoff diagnostic should therefore be interpreted relative to
`q_*`, not only relative to a fixed number in units of `k_F`.
