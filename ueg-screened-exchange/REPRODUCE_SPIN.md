# Reproducing the spin paper

The manuscript source is distributed separately as `prb_manuscript_source.tar.gz`. Its integral-based
figures (3, 4, 5, 7 and 8), four tables, and associated numerical statements
use `production/raw/`. SM Figure S5 is a covariance-aware reanalysis of the saved
spin-map replicates. SM Figure S4 uses the analytic asymptotic coefficients;
Figures 1 and 3 are manuscript-side diagrams; Figure 2 is an analytic profile; Figure 6 is produced from the high-density validation data (see `validation/README.md`). The separate figure-generation
steps below distinguish these data sources.

The commands below are separate steps so that each input and output can be
inspected. Run them from the directory containing this file.

The current notation uses `nu` for imaginary frequency, `L` for the Lindhard
shape, `q_tf_sq` for the squared reference Thomas--Fermi wavevector, and
`W_lambda = lambda*v/epsilon_lambda` for the screened interaction. The energy
is `e_kite = 2 K[Wbar]`, where `Wbar = integral_0^1 W_lambda dlambda = v R(c)/2`.
The numerical energies remain in Ry per electron. Archived raw records and
their source hashes retain their historical field names; the current reader
performs the documented compatibility conversion. Renaming notation does not
replace historical measurements with new calculations.

## 1. Environment and inputs

Use Python **3.12.14** and install the pinned packages:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-spin.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

The released run used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and Matplotlib
3.10.8, with threadpoolctl 3.6.0, on Linux x86_64. The interpreter version is
compared exactly, so a different patch release such as 3.12.3 is refused before
any result is written; install 3.12.14 with pyenv, uv or a distribution image.
The same contract gates `spin/run_inprocess.py`, because that notebook rebuilds
the released products. The command-line entry points
check these exact versions and force one numerical-library thread, a fixed
Python hash seed, and fixed output timestamps before computation.
`production/runtime_lock.json` defines this contract. Each raw record stores it and source
hashes. All integral inputs are in `production/config.json`; the two radial
grids are `data/appendixB_grid_published.dat` and
`data/appendixB_grid_extended.dat`. `production/reference_brm.json` contains
the externally published comparison values; they are not computed by this code.

| Dataset | Figure/table | Samples per term or node, per replicate | Replicates | Other settings |
|---|---|---:|---:|---|
| `tableI` | Fig. 5, SM Fig. S2; Table I | 2^20 | 64 | q scale q*, 6 densities |
| `tableII` | Table III | 2^19 | 8 | q scale 2, 10 densities |
| `tableIII_direct` | SM Table SXVI, direct | 2^21 | 12 | q scale 2, 3 densities |
| `tableIV` | SM Table SXX | 2^19 | 16 | q scale 2, 5 densities |
| `map` | Table IV, SM Figs. S3 and S5 | 2^17 | 8 | q scale 2, 13×8 nodes |
| `shell` | Fig. 4 | 2^15 | 8 | 3×16 nodes; frequency orders 24 and 36 |
| `grid` | SM Table SXVI, quadrature | 15000 | 8 | 65 nodes containing the 49-node grid; frequency orders 24 and 48 |

Direct datasets use master seed 20260713. The shell uses 714203; grid uses 41.
The exact replicate/node seed rules are explicit in `production/compute.py`.
Direct sampling processes 2^15 points at a time to bound memory use while
retaining the full Sobol sequence. The two direct terms use independent
scrambles. Map nodes share the same set of scrambles, preserving covariance.

## 2. Check or redraw the released measurements

This step uses saved replicate means. It does not evaluate an integral.

```bash
python production/verify.py --manifest
python production/verify.py --render
```

The second command generates the original production outputs in a temporary
directory and compares them with the release. It checks their generated PDFs,
TeX tables, numerical-text macros, CSVs, JSON summaries and compatibility data
files. SM Figures S4 and S5 have the separate generation and verification steps below. Under the recorded
plotting environment the files agree byte for byte; PDF timestamps are fixed.

To regenerate the displayed products from the validated released measurements:

```bash
python production/render.py --output ../prb_regenerated_products
```

This writes the `spin/figs/`, `spin/generated/`, `production/generated/`, and
compatibility products beneath `../prb_regenerated_products/`. Do not edit those outputs
by hand. Every staged product is checked against `production/product_lock.json`
before any destination is replaced; a changed figure or table stops the render.
The historical `figs/make_shell_fig.py`, `make_map_figs.py`, and
`make_companion_figs.py` entry points now call this generator.

## 3. Recompute the measurements from their defining integrals

Each command below evaluates one dataset. `--force` means that the integral is
actually rerun, even when a current record exists. Without it, an existing
record with exact or explicitly registered compatible source hashes and the
reference digest is reported as `CURRENT` and skipped; unknown or altered
records cause an error. At the default output location, `--force` checks a
fresh evaluation and preserves the locked historical record (`VERIFIED`).
Use `--output` to retain a new record with current source hashes.

```bash
python production/compute.py --dataset tableI --force
python production/compute.py --dataset tableII --force
python production/compute.py --dataset tableIII_direct --force
python production/compute.py --dataset tableIV --force
python production/compute.py --dataset map --force
python production/compute.py --dataset shell --force
python production/compute.py --dataset grid --force
```

For a selected density, append e.g. `--rs 1000`. For a selected map node, append
`--rs 380 --zeta 1`. A new output directory protects the released measurements:

```bash
python production/compute.py --dataset tableII --rs 5 --output rerun_raw
python production/verify.py --recompute tableII --rs 5
python production/verify.py --recompute shell --rs 1000
python production/verify.py --recompute grid --rs 1
python production/verify.py --recompute map --rs 380 --zeta 1
```

The verification commands recompute in temporary directories and compare
**every numerical value, setting, replicate, and non-source metadata field
exactly**. The historical source hashes remain in the locked records; a
fresh result records the hashes of the code actually used. The only permitted
source differences are the complete old/current hash maps registered in
`production/source_compatibility.json`. Unknown revisions are rejected.
Each fresh record is checked against `production/reference_lock.json` before
being saved to a separate output directory. A mismatch stops execution. Seeds are
attached to terms/replicates/nodes, not worker scheduling; changing the order
of independent jobs cannot change their sequences. Integer and decimal CLI
spellings (e.g. `5` and `5.0`) use the same canonical input representation.
Wall-clock timings go to the console only, never into result JSON.

For a complete independent two-pass check (134 records in each pass):

```bash
python production/repeat_all.py --jobs 2 --passes 2 --report repeat_report.json
```

This evaluates the actual integrals twice in empty temporary directories;
it does not satisfy the check by loading cached measurements. `--jobs` controls
independent processes, each with one BLAS thread. Every pass must also match
the released reference under the same explicit source-compatibility rule,
and the fresh passes must also agree with each other byte for byte.
The report contains hashes, not varying timestamps.

The tested contract is exact results on the recorded runtime. A different CPU,
BLAS build, operating system, or plotting backend can still change floating-point
or PDF bytes even with the same package versions. Such changes are detected,
not silently accepted. This is not a promise of identical arithmetic on every
possible machine. Preserve the released lock and investigate any mismatch.

The strong-screening table is the largest calculation (about 805 million
sample evaluations over its six densities). In the release run each such
density took about four minutes; ordinary-density points took about 13 seconds.
Times depend on hardware. The dataset commands are independent and can be
scheduled separately, with memory use controlled by the chunk size.

Raw files contain independently randomized replicate means, not every
individual sample. The seed, generator, settings and source hashes specify how
to recreate the sample sequence. The raw files for the shell and radial grid
also contain node-by-node replicate values at both frequency orders.

## 4. Export figures and tables; manuscript compilation is separate

```bash
python production/render.py --output ../prb_regenerated_products
python production/verify.py --render
```

Rendering writes the figure PDFs and generated TeX table/numerical-text exports.
Generated TeX tables and numerical macros included in this code archive are
numeric inputs, not manuscript source; the article and derivation-note TeX
remain in the separate manuscript archive. `verify.py --render` regenerates
the original production exports in a temporary directory and checks their
reference digests even when those saved exports are absent. Other missing
generated outputs remain errors. The additional SM Figure S5 numerical macros are
checked by `spin/verify_spin_bias.py`, as described below.

Compile the paper using the instructions in `prb_manuscript_source.tar.gz`.
Copy newly exported `../prb_regenerated_products/spin/figs/` and
`../prb_regenerated_products/spin/generated/` into the manuscript directory
if you intentionally regenerate them. Numerical computation and the notebook
require no TeX installation.

## 5. Execute the reviewer notebook

```bash
cd spin
python run_inprocess.py
python -m nbconvert --to html spin_reviewer_verification.ipynb
cd ..
```

Alternatively open the notebook in Jupyter, restart the kernel and Run All.
Run it in `spin/` with this code/data bundle present. Manuscript-source checks
are explicitly skipped because this distribution excludes the manuscript. Default execution includes
the independent imaginary-axis integrator and verifies saved production
records and regenerated products. It does not repeat all production integrals;
`RUN_PRODUCTION=True` explicitly enables that additional calculation.

The notebook separates symbolic identities, numerical quadrature, model
assembly, data transformations, and fresh original-integral diagnostics.
Reading a stored measurement is not counted as an independent integral test.

## 6. Changes and future releases

Normal reruns never update the reference lock. If an intentional algorithm or
configuration change requires a new scientific baseline, first write candidates
outside the released raw directory, explicitly using the maintainer option:

```bash
python production/compute.py --dataset tableII --prepare-release --output candidate_raw
```

Repeat this for **every dataset** because source provenance is part of every
record. Review numerical changes, repeat all candidate calculations independently,
and only then promote the candidate records and their SHA-256 values to
`production/reference_lock.json`. Review and update `production/product_lock.json`
for any intentionally changed figures/tables as a separate release action. This is a manual scientific release decision;
there is no automatic rebaseline command. `--prepare-release` is rejected when
its output is the released `production/raw` directory.

After promotion, render the outputs, compile the manuscript, execute the
notebook, and review the figures and numerical changes. Then seal integrity:

```bash
python production/verify.py --render
python production/verify.py --seal
python production/verify.py --manifest
```

The manifest records the integrity of the reviewed release; it does not prove
the mathematical derivation. `--seal` is an explicit release step, not an
automatic repair operation. The checks in `kernels_checks/` and `kernels_review/` are documented separately;
this document governs the spin paper in `spin/`.

Guard regression checks (no expensive integral):

```bash
python production/test_determinism.py
```

## Analytic spin-interpolation figure (FRG-context revision)

SM Figure S4 plots the ratio of the three-term analytic energies and its first
correction. It uses no Monte Carlo data and makes no FRG prediction. Run:

```bash
python figs/make_spin_interpolation.py
python spin/verify_spin_interpolation.py
```

The PDF is `figs/fig_spin_interpolation.pdf`; plotted values are
`figs/data/fig_spin_interpolation.csv`. Copy the PDF into the manuscript's
`figs/` directory after regeneration. Check 32 of the reviewer notebook
executes the same verification. The original numerical payloads are preserved. Source-path keys and their file locks were normalized as documented in the current migration receipt. The matching schematic remains Figure 2.

## Spin-interpolation bias from saved integral data (SM Figure S5)

Run these steps from the code root:

```bash
python figs/make_spin_bias.py
python spin/verify_spin_bias.py
python spin/verify_spin_bias_algebra.py
```

The generator reads all 104 saved `map` records and forms the endpoint
interpolation residual separately for each of the eight shared scrambles and
each of the two independently sampled terms. The standard error retains the
covariance among the intermediate polarization and both endpoints; variances
of the two independent terms are then added. This step performs no new
Monte Carlo integral and fits no asymptotic coefficient.

Outputs are `figs/fig_spin_bias.pdf`, `figs/data/fig_spin_bias.csv` and
`spin/generated/spin_bias_numbers.tex`. Copy the PDF and macro input into the
corresponding `figs/` and `generated/` directories of the separate manuscript
if you intentionally regenerate them. The CSV retains all 104 residual rows;
the figure displays two intermediate polarizations. The macro file contains
only numeric inputs for the article's three examples, not manuscript prose.

The first check independently propagates the full three-polarization
covariance matrices and verifies the exported numbers and reproducible
outputs. The second checks the interpolation algebra, the analytic sign proof,
and 70-digit supplementary evaluations. Neither replaces a fresh integral
calculation or proves statistical confidence coverage. The current derivation
notes cover all 88 numbered displays.
