#!/usr/bin/env bash
# Rebuild tables and figures from the released, already-computed raw results.
# This is deliberately separate from reproduce.sh, which executes ks.x.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
OUTPUT_ROOT=""

usage() {
  cat <<'EOF'
Usage: ./postprocess_results.sh [--output DIR]

Regenerate Tables II/III, validate their rounded manuscript values, regenerate
all six figures, and parse the constant-kF benchmark from committed data under
results/raw.  No ks.x calculation is performed.  Outputs are written to an
untracked directory so the released files remain unchanged.

Options:
  --output DIR   Output directory.  Default: runs/postprocess/<UTC timestamp>/
  -h, --help     Show this help.
EOF
}

while (($#)); do
  case "$1" in
    --output)
      (($# >= 2)) || { echo "ERROR: --output requires a directory" >&2; exit 2; }
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option '$1'" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "$OUTPUT_ROOT" ]]; then
  OUTPUT_ROOT="$ROOT/runs/postprocess/$(date -u +%Y%m%dT%H%M%SZ)"
elif [[ "$OUTPUT_ROOT" != /* ]]; then
  OUTPUT_ROOT="$ROOT/$OUTPUT_ROOT"
fi
command -v realpath >/dev/null 2>&1 || {
  echo "ERROR: realpath is required for safe output-path validation" >&2
  exit 2
}
OUTPUT_ROOT="$(realpath -m -- "$OUTPUT_ROOT")"

# Inside this repository, generated files belong only below runs/.  Absolute
# paths outside the repository (for example /tmp/postprocess) are also allowed.
case "$OUTPUT_ROOT/" in
  "$ROOT/runs/"*) ;;
  "$ROOT/"*)
    echo "ERROR: repository-internal output must be below $ROOT/runs" >&2
    exit 2
    ;;
esac

if [[ -d "$OUTPUT_ROOT" ]] && find "$OUTPUT_ROOT" -mindepth 1 -print -quit | grep -q .; then
  echo "ERROR: refusing to overwrite nonempty directory: $OUTPUT_ROOT" >&2
  exit 2
fi
mkdir -p "$OUTPUT_ROOT/tables" "$OUTPUT_ROOT/figures"
OUTPUT_ROOT="$(cd "$OUTPUT_ROOT" && pwd)"

echo '[tables] rebuilding from committed results/raw'
python3 "$ROOT/scripts/make_tables.py" \
  --raw "$ROOT/results/raw" --ref "$ROOT/results/reference" \
  --out "$OUTPUT_ROOT/tables"

echo '[tables] validating printed manuscript values'
python3 "$ROOT/scripts/validate_tables.py" \
  --tables "$OUTPUT_ROOT/tables" \
  --expected "$ROOT/results/reference/paper_example_values.json"

echo '[figures] regenerating all six figures from committed data'
python3 "$ROOT/scripts/make_release_figures.py" \
  --gp-dir "$ROOT/results/gp" --output-dir "$OUTPUT_ROOT/figures"

echo '[benchmark] parsing committed constant-kF diagonalization results'
python3 "$ROOT/scripts/make_diag_benchmark.py" \
  --raw "$ROOT/results/raw/diag_benchmark" \
  | tee "$OUTPUT_ROOT/diag_benchmark.txt"

echo "DONE: postprocessed artifacts are in $OUTPUT_ROOT"
