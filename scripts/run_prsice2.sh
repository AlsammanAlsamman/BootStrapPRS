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
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -n "$prsice_bin" && -n "$prsice_r" && -n "$rscript_cmd" && -n "$base" && -n "$target" && -n "$out" && -n "$done" ]] || usage

[[ -f "$prsice_bin" ]] || { echo "Missing PRSice binary: $prsice_bin" >&2; exit 1; }
[[ -f "$prsice_r" ]] || { echo "Missing PRSice.R: $prsice_r" >&2; exit 1; }
[[ -f "$base" ]] || { echo "Missing base: $base" >&2; exit 1; }
[[ -f "${target}.bed" && -f "${target}.bim" && -f "${target}.fam" ]] || { echo "Missing target PLINK files for prefix: $target" >&2; exit 1; }

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

echo "[run_prsice2] ${rscript_cmd} ${prsice_r} ${args[*]}" >&2
"$rscript_cmd" "$prsice_r" "${args[@]}"

mkdir -p "$(dirname "$done")"
echo "ok" > "$done"
