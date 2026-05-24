#!/usr/bin/env bash
set -euo pipefail

# Make the build independent of the caller's current working directory and
# keep all compiler products inside src/ (which is ignored by Git).
SRC_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SRC_DIR"

gfortran -O2 -ffree-line-length-none \
  nofdft_ks_spherical_skp_nlf_v8.f90 -o ks.x -llapack -lblas
echo "built $SRC_DIR/ks.x"
