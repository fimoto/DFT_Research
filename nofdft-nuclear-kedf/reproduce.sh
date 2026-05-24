#!/usr/bin/env bash
# Prepare or run the released self-consistent OFDFT calculation(s).  Archived
# results are not acceptance targets: comparison with them is an explicit,
# optional diagnostic.  Repository-internal output is confined to runs/.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd -P)"
TAGS=(O16 Ca40 A80 A140)
CASES=()
RUN_ALL=0
BUILD=1
COMPARE_ARCHIVED=0
PREPARE_ONLY=0
OUTPUT_ROOT=""

usage() {
  cat <<'EOF'
Usage:
  ./reproduce.sh --case FAMILY/NUCLEUS [--case ...] [options]
  ./reproduce.sh --all [options]
  ./reproduce.sh --list

Run the released ks.x solver from an exact production namelist, starting from
the released KS species density.  Each selected case performs three stages:
an nmu=24 self-consistent ITE relaxation, an nmu=64 energy reevaluation, and an
nmu=48 NLF reevaluation.  New outputs are collected in an isolated run
directory and are not compared with archived results by default.

Families: per-nucleus, shared
Nuclei:   O16, Ca40, A80, A140

Options:
  --case CASE       Run one case, e.g. per-nucleus/O16 (repeatable).
  --all             Run all eight family/nucleus cases.
  --output DIR      Parent output directory.  Default:
                    runs/reproduction/<UTC timestamp>/
  --no-build        Use the existing src/ks.x instead of rebuilding it.
  --prepare-only    Stage all exact inputs and provenance records, but do not
                    build or execute ks.x.
  --compare-archived
                    After a successful run, compare with the historical files
                    in results/raw/ (optional diagnostic, not an acceptance
                    test for newly generated production data).
  --list            List valid cases and exit.
  -h, --help        Show this help.

Relative --output paths are resolved from the repository root and must lie
under runs/.  Absolute paths outside the repository are also allowed.  Existing
nonempty case directories are never overwritten.
EOF
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

list_cases() {
  local family tag
  for family in per-nucleus shared; do
    for tag in "${TAGS[@]}"; do
      printf '%s/%s\n' "$family" "$tag"
    done
  done
}

normalize_case() {
  local value="$1" family tag
  family="${value%%/*}"
  tag="${value#*/}"
  [[ "$value" == */* && -n "$family" && -n "$tag" ]] || return 1
  case "$family" in
    per-nucleus|per_nucleus) family="per-nucleus" ;;
    shared) ;;
    *) return 1 ;;
  esac
  case "$tag" in
    O16|Ca40|A80|A140) ;;
    *) return 1 ;;
  esac
  printf '%s/%s\n' "$family" "$tag"
}

while (($#)); do
  case "$1" in
    --case)
      (($# >= 2)) || die "--case requires FAMILY/NUCLEUS"
      normalized="$(normalize_case "$2")" || die "invalid case '$2' (use --list)"
      CASES+=("$normalized")
      shift 2
      ;;
    --all) RUN_ALL=1; shift ;;
    --output)
      (($# >= 2)) || die "--output requires a directory"
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --no-build) BUILD=0; shift ;;
    --prepare-only) PREPARE_ONLY=1; shift ;;
    --compare-archived) COMPARE_ARCHIVED=1; shift ;;
    --list) list_cases; exit 0 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option '$1' (use --help)" ;;
  esac
done

if ((RUN_ALL)); then
  ((${#CASES[@]} == 0)) || die "use either --all or --case, not both"
  for family in per-nucleus shared; do
    for tag in "${TAGS[@]}"; do
      CASES+=("$family/$tag")
    done
  done
fi
((${#CASES[@]} > 0)) || { usage >&2; exit 2; }

if ((PREPARE_ONLY && COMPARE_ARCHIVED)); then
  die "--prepare-only cannot be combined with --compare-archived"
fi
if ((COMPARE_ARCHIVED)); then
  command -v python3 >/dev/null 2>&1 || die "python3 is required for result comparison"
fi
command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required for provenance records"
command -v realpath >/dev/null 2>&1 || die "realpath is required for safe output-path handling"

if [[ -z "$OUTPUT_ROOT" ]]; then
  OUTPUT_ROOT="$ROOT/runs/reproduction/$(date -u +%Y%m%dT%H%M%SZ)"
elif [[ "$OUTPUT_ROOT" != /* ]]; then
  OUTPUT_ROOT="$ROOT/$OUTPUT_ROOT"
fi
OUTPUT_ROOT="$(realpath -m -- "$OUTPUT_ROOT")"

# A relative path is repository-internal, and an absolute path may resolve
# there through symlinks.  In either case, confine repository-internal writes
# to the ignored runs/ subtree.  This protects source, inputs, and results/raw.
case "$OUTPUT_ROOT/" in
  "$ROOT/"|"$ROOT/"*)
    case "$OUTPUT_ROOT/" in
      "$ROOT/runs/"*) ;;
      *) die "repository-internal --output must be under $ROOT/runs" ;;
    esac
    ;;
esac
mkdir -p "$OUTPUT_ROOT"

RUN_INDEX="$OUTPUT_ROOT/run_index.tsv"
if [[ ! -e "$RUN_INDEX" ]]; then
  printf 'case_id\tstatus\ttimestamp_utc\telapsed_seconds\tcase_directory\n' >"$RUN_INDEX"
elif [[ ! -f "$RUN_INDEX" ]]; then
  die "run index path is not a regular file: $RUN_INDEX"
fi

if ((!PREPARE_ONLY && BUILD)); then
  printf '[build] %s\n' "$ROOT/src/ks.x"
  (cd "$ROOT/src" && ./build.sh)
fi
BIN="$ROOT/src/ks.x"
if ((!PREPARE_ONLY)); then
  [[ -x "$BIN" ]] || die "$BIN is missing or not executable; build it or omit --no-build"
fi

record_index() {
  local requested="$1" status="$2" timestamp="$3" elapsed="$4" case_dir="$5"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$requested" "$status" "$timestamp" "$elapsed" "$case_dir" >>"$RUN_INDEX"
}

write_metadata() {
  local file="$1" requested="$2" status="$3" prepared="$4"
  local started="$5" finished="$6" elapsed="$7" density_source="$8" kernel_source="$9"
  {
    printf 'case=%s\n' "$requested"
    printf 'status=%s\n' "$status"
    printf 'prepared_utc=%s\n' "$prepared"
    [[ -z "$started" ]] || printf 'started_utc=%s\n' "$started"
    [[ -z "$finished" ]] || printf 'finished_utc=%s\n' "$finished"
    [[ -z "$elapsed" ]] || printf 'elapsed_seconds=%s\n' "$elapsed"
    printf 'solver=src/ks.x\n'
    printf 'source_density=%s\n' "${density_source#"$ROOT/"}"
    printf 'source_kernel=%s\n' "${kernel_source#"$ROOT/"}"
    printf 'system=%s\n' "$(uname -srm)"
    if command -v gfortran >/dev/null 2>&1; then
      printf 'compiler=%s\n' "$(gfortran --version | head -n 1)"
    fi
  } >"$file"
}

stage_error() {
  local case_dir="$1" stage="$2" message="$3"
  printf '  [%s] INVALID OUTPUT: %s\n' "$stage" "$message" >&2
  printf '  [%s] last stdout/stderr lines follow\n' "$stage" >&2
  tail -n 20 "$case_dir/$stage.stdout.log" >&2 || true
  tail -n 20 "$case_dir/$stage.stderr.log" >&2 || true
  return 1
}

validate_stage() {
  local case_dir="$1" stage="$2" required=() file
  if ! grep -Fq 'All requested tasks finished.' "$case_dir/$stage.stdout.log"; then
    stage_error "$case_dir" "$stage" "solver completion sentinel is absent"
    return 1
  fi
  if grep -Fq 'ERROR:' "$case_dir/$stage.stdout.log" || \
     grep -Fq 'ERROR:' "$case_dir/$stage.stderr.log"; then
    stage_error "$case_dir" "$stage" "solver log contains ERROR:"
    return 1
  fi

  case "$stage" in
    ite)
      required=(ite_history.dat ite_density.dat ite_derivative_radial.dat ite_summary.dat ite_kernel_radial.dat)
      ;;
    energy)
      required=(energy.dat derivative_radial.dat derivative_kernel.dat)
      ;;
    nlf)
      required=(nlf.dat)
      ;;
    *) die "internal error: unknown stage '$stage'" ;;
  esac
  for file in "${required[@]}"; do
    if [[ ! -s "$case_dir/$file" ]]; then
      stage_error "$case_dir" "$stage" "missing or empty $file"
      return 1
    fi
  done

  case "$stage" in
    ite)
      if ! grep -Eq '^#[[:space:]]+attempted_iterations[[:space:]]+8000[[:space:]]*$' \
          "$case_dir/ite_history.dat"; then
        stage_error "$case_dir" "$stage" \
          "ITE history lacks the final '# attempted_iterations 8000' footer"
        return 1
      fi
      if ! grep -Fq 'E_total_NeqZ_MeV' "$case_dir/ite_summary.dat"; then
        stage_error "$case_dir" "$stage" "ite_summary.dat lacks E_total_NeqZ_MeV"
        return 1
      fi
      if ! awk '
        /^[[:space:]]*#/ || NF == 0 { next }
        { rows++; if (NF < 2) bad = 1 }
        END { exit ! (rows >= 250 && !bad) }
      ' "$case_dir/ite_density.dat"; then
        stage_error "$case_dir" "$stage" \
          "ite_density.dat must contain at least 250 data rows with at least 2 columns"
        return 1
      fi
      ;;
    energy)
      if ! grep -Fq 'E_total_NeqZ_MeV' "$case_dir/energy.dat"; then
        stage_error "$case_dir" "$stage" "energy.dat lacks E_total_NeqZ_MeV"
        return 1
      fi
      ;;
    nlf)
      if ! awk '
        /^[[:space:]]*#/ || NF == 0 { next }
        { rows++; if (NF < 7) bad = 1 }
        END { exit ! (rows >= 250 && !bad) }
      ' "$case_dir/nlf.dat"; then
        stage_error "$case_dir" "$stage" \
          "nlf.dat must contain at least 250 data rows with at least 7 columns"
        return 1
      fi
      ;;
  esac
}

run_stage() {
  local case_dir="$1" stage="$2"
  cp "$case_dir/$stage.nml" "$case_dir/input_ks_nlf.dat"
  printf '  [%s] running ks.x\n' "$stage"
  if ! (cd "$case_dir" && "$BIN" >"$stage.stdout.log" 2>"$stage.stderr.log"); then
    printf '  [%s] FAILED; last stdout/stderr lines follow\n' "$stage" >&2
    tail -n 20 "$case_dir/$stage.stdout.log" >&2 || true
    tail -n 20 "$case_dir/$stage.stderr.log" >&2 || true
    return 1
  fi
  validate_stage "$case_dir" "$stage"
}

run_case() {
  local requested="$1" family="${1%%/*}" tag="${1#*/}"
  local family_dir kernel_source density_source case_dir input_dir prepared started finished
  local start_s end_s elapsed status
  if [[ "$family" == "per-nucleus" ]]; then
    family_dir="per_nucleus"
    kernel_source="$ROOT/kernels/per_nucleus/eg_PN14_${tag}.dat"
  else
    family_dir="shared"
    kernel_source="$ROOT/kernels/universal/eg_SH14.dat"
  fi
  input_dir="$ROOT/production_inputs/$family_dir/$tag"
  case_dir="$OUTPUT_ROOT/$family_dir/$tag"
  density_source="$ROOT/inputs/ks/ks_${tag}_nr250.dat"

  [[ -d "$input_dir" ]] || die "missing production inputs: $input_dir"
  [[ -f "$density_source" ]] || die "missing KS density for $tag"
  [[ -f "$kernel_source" ]] || die "missing kernel: $kernel_source"
  if [[ -d "$case_dir" ]] && find "$case_dir" -mindepth 1 -print -quit | grep -q .; then
    die "refusing to overwrite nonempty case directory: $case_dir"
  fi
  mkdir -p "$case_dir"

  cp "$density_source" "$case_dir/initial_density.dat"
  cp "$kernel_source" "$case_dir/kernel.dat"
  cp "$input_dir/ite.nml" "$input_dir/energy.nml" "$input_dir/nlf.nml" "$case_dir/"

  (
    cd "$case_dir"
    sha256sum initial_density.dat kernel.dat ite.nml energy.nml nlf.nml
  ) >"$case_dir/input_sha256.txt"

  prepared="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  write_metadata "$case_dir/run_metadata.txt" "$requested" prepared "$prepared" "" "" "" \
    "$density_source" "$kernel_source"
  record_index "$requested" prepared "$prepared" "" "$family_dir/$tag"
  printf '[case] %s -> %s\n' "$requested" "$case_dir"

  if ((PREPARE_ONLY)); then
    printf '  [prepared] exact inputs staged; solver not executed\n'
    return
  fi

  (cd "$ROOT" && sha256sum src/ks.x) >"$case_dir/solver_sha256.txt"
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  start_s="$(date +%s)"
  if run_stage "$case_dir" ite && \
     run_stage "$case_dir" energy && \
     run_stage "$case_dir" nlf; then
    status="completed"
  else
    status="failed"
  fi
  end_s="$(date +%s)"
  elapsed="$((end_s - start_s))"
  finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  write_metadata "$case_dir/run_metadata.txt" "$requested" "$status" "$prepared" "$started" \
    "$finished" "$elapsed" "$density_source" "$kernel_source"
  record_index "$requested" "$status" "$finished" "$elapsed" "$family_dir/$tag"
  [[ "$status" == completed ]] || return 1

  if ((COMPARE_ARCHIVED)); then
    printf '  [compare] historical archived output (optional diagnostic)\n'
    python3 "$ROOT/scripts/compare_reproduction.py" \
      --case "$requested" --run "$case_dir" \
      --repository-root "$ROOT" --json "$case_dir/comparison.json" \
      | tee "$case_dir/comparison.txt"
  fi
}

for item in "${CASES[@]}"; do
  run_case "$item"
done

if ((PREPARE_ONLY)); then
  printf 'DONE: exact calculation inputs are staged in %s\n' "$OUTPUT_ROOT"
else
  printf 'DONE: solver outputs and logs are in %s\n' "$OUTPUT_ROOT"
fi
