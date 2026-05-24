# Ten-density comparison of the independent implementation

`independent_audit.ipynb` re-evaluates the kite energy at three densities, which
is what fits in the notebook's six-minute budget. The ten-density comparison
quoted in the Supplemental Material was produced by the scripts here, with the
same estimator and the same module (`a7_mc.py`), at larger sample sizes.

| script | what it produced | recorded output |
|---|---|---|
| `brm_compare.py` | r_s = 0.1 to 30 against spin-paper Table II and BRM Table IX | `brm_cmp.log` |
| `rs50.py` | the r_s = 50 row, which the first run did not reach | `rs50.log` |
| `ref30.py` | the r_s = 30 grid-refinement pair, 84 and 160 momentum nodes | `ref30.log` |
| `table_S_grid.py` | SM Table S-grid, r_s = 1e8 under three momentum grids, plus two convergence checks | `table_S_grid.json` (every replicate at full precision), `table_S_grid.log` |

Settings are in the scripts, not in prose: eight replicates (`brm_compare.py`,
`rs50.py`) or six (`ref30.py`, whose standard error uses sqrt 6), `N = 200000` or
`250000` points per replicate per term, seeds `3000+i` and `7000+i`, the split
momentum grid of `hd_split.energy_split`, and 40 to 56 frequency nodes. The logs
hold the means and standard errors as printed, and `rs50.log` also lists its
eight replicate values rounded to three decimals; full-precision replicate values
were not retained, so these runs are reproducible from the scripts but not
auditable replicate by replicate at full precision. That is the difference between these numbers and the
spin-paper production records, which do store every replicate.

    python3 brm_compare.py     # about 40 minutes on one core
    python3 rs50.py
    python3 ref30.py

`a7_mc.py` and `hd_split.py` are the estimator; they import nothing from
`src/`, `production/` or `pysecdec_kite/`.

`table_S_grid.py` (added 2026-09-16) replaces values whose script and output were
not in the archive. It uses eight replicates of 10^6 samples per momentum node and
stores every replicate value, so unlike the runs above it is auditable replicate
by replicate. At 2.5e5 samples the replicate distribution at r_s = 1e8 has a heavy
tail and the eight-replicate error understated the scatter; the checks in the same
file (a finer momentum split and frequency order 80) agree with the table rows and
with the three-term expansion within one standard error.

    python3 table_S_grid.py     # about 40 minutes on one core
