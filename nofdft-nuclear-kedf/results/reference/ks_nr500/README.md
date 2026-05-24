# ks_nr500: raw outputs of the Kohn-Sham reference calculations (STEP 1)

Species densities and single-particle levels of the four
nr=500 KS reference runs (namelists/ks/ks_<TAG>_nr500.nml). These are the
regenerated outputs of the released namelists, archived for completeness; the
production chain uses the FROZEN derived products from the original runs:
inputs/ks/ks_<TAG>_nr250.dat (resampled input densities) and
results/reference/ks_nlf_<TAG>.dat (tau/NLF reference tables). Agreement of
the regenerated solutions with the frozen products: max|drho| ~ 6e-10 for
O16/Ca40 and ~3e-6 fm^-3 for A80/A140 (SCF tolerance 1e-8; independent
convergence paths), i.e. far below every reported observable's precision.
