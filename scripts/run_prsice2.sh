#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Run PRSice-2 using a formatted base file and a PLINK target prefix.

Required:
  --prsice-bin PATH
  --prsice-r   PATH
  --rscript    CMD
  --base       FILE
  --target     PREFIX
  --out        PREFIX
  --done       FILE

Optional:
  --ld PREFIX
  --stat OR|BETA
  --score avg|sum|std|con-std
  --binary-target T|F
  --clump-kb N
  --clump-r2 X
  --clump-p  X
  --bar-levels "..."
  --perm N
  --thread N
  --quantile N
  --quant-break "..."
  --quant-ref N
  --extra "<raw extra args>"
EOF
  exit 2
}

prsice_bin=""
prsice_r=""
rscript_cmd=""
base=""
target=""
out=""
done=""

ld=""
stat="OR"
score="avg"
binary_target=""
clump_kb=""
clump_r2=""
clump_p=""
bar_levels=""
perm=""
thread=""
quantile=""
quant_break=""
quant_ref=""
extra=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prsice-bin) prsice_bin="$2"; shift 2;;
    --prsice-r) prsice_r="$2"; shift 2;;
    --rscript) rscript_cmd="$2"; shift 2;;
    --base) base="$2"; shift 2;;
    --target) target="$2"; shift 2;;
    --out) out="$2"; shift 2;;
    --done) done="$2"; shift 2;;
    --ld) ld="$2"; shift 2;;
    --stat) stat="$2"; shift 2;;
    --score) score="$2"; shift 2;;
    --binary-target) binary_target="$2"; shift 2;;
    --clump-kb) clump_kb="$2"; shift 2;;
    --clump-r2) clump_r2="$2"; shift 2;;
    --clump-p) clump_p="$2"; shift 2;;
    --bar-levels) bar_levels="$2"; shift 2;;
    --perm) perm="$2"; shift 2;;
    --thread) thread="$2"; shift 2;;
    --quantile) quantile="$2"; shift 2;;
    --quant-break) quant_break="$2"; shift 2;;
    --quant-ref) quant_ref="$2"; shift 2;;
    --extra) extra="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -n "$prsice_bin" && -n "$prsice_r" && -n "$rscript_cmd" && -n "$base" && -n "$target" && -n "$out" && -n "$done" ]] || usage

[[ -f "$prsice_bin" ]] || { echo "Missing PRSice binary: $prsice_bin" >&2; exit 1; }
[[ -f "$prsice_r" ]] || { echo "Missing PRSice.R: $prsice_r" >&2; exit 1; }
[[ -f "$base" ]] || { echo "Missing base: $base" >&2; exit 1; }

# Best-effort: make PRSice binary executable on POSIX filesystems
chmod +x "$prsice_bin" 2>/dev/null || true

if ! command -v "$rscript_cmd" >/dev/null 2>&1; then
  echo "Rscript command not found in PATH: $rscript_cmd" >&2
  echo "Hint: set prsice2.module in analysis.yml to load R" >&2
  exit 1
fi

if [[ ! -f "${target}.bed" || ! -f "${target}.bim" || ! -f "${target}.fam" ]]; then
  echo "Missing target PLINK files for prefix: $target" >&2
  ls -l "${target}.bed" "${target}.bim" "${target}.fam" 2>/dev/null || true
  exit 1
fi

mkdir -p "$(dirname "$out")"

args=(
  --prsice "$prsice_bin"
  --base "$base"
  --target "$target"
  --out "$out"
  --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --stat "$stat" --pvalue P
  --score "$score"
)

if [[ -n "$binary_target" && "$binary_target" != "NONE" && "$binary_target" != "None" ]]; then
  args+=(--binary-target "$binary_target")
fi

if [[ -n "$ld" && "$ld" != "NONE" && "$ld" != "None" ]]; then
  args+=(--ld "$ld")
fi

if [[ -n "$clump_kb" && "$clump_kb" != "NONE" && "$clump_kb" != "None" ]]; then
  args+=(--clump-kb "$clump_kb")
fi
if [[ -n "$clump_r2" && "$clump_r2" != "NONE" && "$clump_r2" != "None" ]]; then
  args+=(--clump-r2 "$clump_r2")
fi
if [[ -n "$clump_p" && "$clump_p" != "NONE" && "$clump_p" != "None" ]]; then
  args+=(--clump-p "$clump_p")
fi

if [[ -n "$bar_levels" && "$bar_levels" != "NONE" && "$bar_levels" != "None" ]]; then
  args+=(--bar-levels "$bar_levels")
fi

if [[ -n "$perm" && "$perm" != "NONE" && "$perm" != "None" ]]; then
  args+=(--perm "$perm")
fi

if [[ -n "$thread" && "$thread" != "NONE" && "$thread" != "None" ]]; then
  args+=(--thread "$thread")
fi

if [[ -n "$quantile" && "$quantile" != "NONE" && "$quantile" != "None" ]]; then
  args+=(--quantile "$quantile")
fi

if [[ -n "$quant_break" && "$quant_break" != "NONE" && "$quant_break" != "None" ]]; then
  args+=(--quant-break "$quant_break")
fi

if [[ -n "$quant_ref" && "$quant_ref" != "NONE" && "$quant_ref" != "None" ]]; then
  args+=(--quant-ref "$quant_ref")
fi

if [[ -n "$extra" && "$extra" != "NONE" && "$extra" != "None" ]]; then
  # shellcheck disable=SC2206
  extra_arr=( $extra )
  args+=("${extra_arr[@]}")
fi

echo "[run_prsice2] ${rscript_cmd} ${prsice_r} ${args[*]}" >&2

set +e
"$rscript_cmd" "$prsice_r" "${args[@]}"
rc=$?
set -e

# PRSice can fail during plotting with uneven quantiles. If that happens,
# retry once without --quant-break/--quant-ref (still keeps --quantile).
if [[ $rc -ne 0 && -n "$quant_break" && "$quant_break" != "NONE" && "$quant_break" != "None" ]]; then
  echo "[run_prsice2] PRSice failed (exit=$rc). Retrying without --quant-break/--quant-ref" >&2

  args_retry=()
  skip_next=0
  for ((i=0; i<${#args[@]}; i++)); do
    if (( skip_next )); then
      skip_next=0
      continue
    fi
    case "${args[$i]}" in
      --quant-break|--quant-ref)
        skip_next=1
        ;;
      *)
        args_retry+=("${args[$i]}")
        ;;
    esac
  done

  echo "[run_prsice2] ${rscript_cmd} ${prsice_r} ${args_retry[*]}" >&2
  set +e
  "$rscript_cmd" "$prsice_r" "${args_retry[@]}"
  rc=$?
  set -e
fi

# Final fallback: if plotting still fails, retry without quantile plotting entirely.
if [[ $rc -ne 0 && -n "$quantile" && "$quantile" != "NONE" && "$quantile" != "None" ]]; then
  echo "[run_prsice2] PRSice still failed (exit=$rc). Retrying without --quantile/--quant-break/--quant-ref" >&2

  args_retry2=()
  skip_next=0
  for ((i=0; i<${#args[@]}; i++)); do
    if (( skip_next )); then
      skip_next=0
      continue
    fi
    case "${args[$i]}" in
      --quantile|--quant-break|--quant-ref)
        skip_next=1
        ;;
      *)
        args_retry2+=("${args[$i]}")
        ;;
    esac
  done

  echo "[run_prsice2] ${rscript_cmd} ${prsice_r} ${args_retry2[*]}" >&2
  set +e
  "$rscript_cmd" "$prsice_r" "${args_retry2[@]}"
  rc=$?
  set -e
fi

if [[ $rc -ne 0 ]]; then
  exit $rc
fi

mkdir -p "$(dirname "$done")"
echo "ok" > "$done"
