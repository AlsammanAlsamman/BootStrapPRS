#!/bin/bash
set -euo pipefail

# Splits a PLINK bed set into train/test for a bootstrap replicate,
# extracts covariate PCs for training samples from a precomputed eigenvec,
# then runs GWAS on the training set.
# This script must NOT read YAML; all parameters must be passed via CLI.

usage() {
  echo "Usage: $0 --plink-cmd CMD --bfile PREFIX --pca-eigenvec FILE --out-dir DIR --rep N --seed-base S --train-frac F --pcs K --done PATH [--maf X] [--geno X] [--mind X] [--hwe X] [--memory ARG] [--threads ARG]" >&2
  exit 2
}

PLINK_CMD=""
BFILE=""
PCA_EIGENVEC=""
OUT_DIR=""
REP=""
SEED_BASE=""
TRAIN_FRAC=""
PCS=""
DONE=""
MAF=""
GENO=""
MIND=""
HWE=""
MEMORY=""
THREADS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --plink-cmd) PLINK_CMD="$2"; shift 2;;
    --bfile) BFILE="$2"; shift 2;;
    --pca-eigenvec) PCA_EIGENVEC="$2"; shift 2;;
    --out-dir) OUT_DIR="$2"; shift 2;;
    --rep) REP="$2"; shift 2;;
    --seed-base) SEED_BASE="$2"; shift 2;;
    --train-frac) TRAIN_FRAC="$2"; shift 2;;
    --pcs) PCS="$2"; shift 2;;
    --done) DONE="$2"; shift 2;;
    --maf) MAF="$2"; shift 2;;
    --geno) GENO="$2"; shift 2;;
    --mind) MIND="$2"; shift 2;;
    --hwe) HWE="$2"; shift 2;;
    --memory) MEMORY="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    *) usage;;
  esac
done

[[ -n "$PLINK_CMD" && -n "$BFILE" && -n "$PCA_EIGENVEC" && -n "$OUT_DIR" && -n "$REP" && -n "$SEED_BASE" && -n "$TRAIN_FRAC" && -n "$PCS" && -n "$DONE" ]] || usage

mkdir -p "$OUT_DIR"

seed=$(( SEED_BASE + REP ))

echo "[split_train_pca_gwas.sh] seed=${seed} rep=${REP} train_frac=${TRAIN_FRAC} pcs=${PCS}"
echo "[split_train_pca_gwas.sh] bfile=${BFILE} out_dir=${OUT_DIR}"

hash_line() {
  # reads stdin, prints hex digest
  if command -v md5sum >/dev/null 2>&1; then
    md5sum | awk '{print $1}'
  elif command -v md5 >/dev/null 2>&1; then
    md5 -q
  elif command -v openssl >/dev/null 2>&1; then
    openssl dgst -md5 | awk '{print $2}'
  else
    echo "ERROR: need md5sum (or md5/openssl) for deterministic split" >&2
    return 127
  fi
}

# Read samples from FAM and deterministically shuffle by hashing (seed, FID, IID).
# Then take the first N_train as training.
keep_train="$OUT_DIR/train.keep"
keep_test="$OUT_DIR/test.keep"

fam_file="${BFILE}.fam"
[[ -f "$fam_file" ]] || { echo "Missing: $fam_file" >&2; exit 1; }

n_total=$(awk 'NF>=2{c++} END{print c+0}' "$fam_file")
if [[ "$n_total" -lt 2 ]]; then
  echo "Not enough samples in $fam_file" >&2
  exit 1
fi

# Round to nearest int, and ensure at least 1 and at most n_total-1
n_train=$(awk -v n="$n_total" -v f="$TRAIN_FRAC" 'BEGIN{v=int(n*f+0.5); if(v<1)v=1; if(v>n-1)v=n-1; print v}')

tmp_sorted="$OUT_DIR/_ids_sorted.txt"
awk '{print $1"\t"$2}' "$fam_file" | \
  while IFS=$'\t' read -r fid iid; do
    h=$(printf '%s\t%s\t%s\n' "$seed" "$fid" "$iid" | hash_line)
    printf '%s\t%s\t%s\n' "$h" "$fid" "$iid"
  done | sort -k1,1 > "$tmp_sorted"

head -n "$n_train" "$tmp_sorted" | awk '{print $2"\t"$3}' > "$keep_train"
tail -n $(( n_total - n_train )) "$tmp_sorted" | awk '{print $2"\t"$3}' > "$keep_test"

# Create separate PLINK bed sets for train/test
train_prefix="$OUT_DIR/train"
test_prefix="$OUT_DIR/test"

$PLINK_CMD --bfile "$BFILE" --keep "$keep_train" --make-bed --allow-no-sex ${MEMORY:-} ${THREADS:-} --out "$train_prefix"
$PLINK_CMD --bfile "$BFILE" --keep "$keep_test" --make-bed --allow-no-sex ${MEMORY:-} ${THREADS:-} --out "$test_prefix"

# Build covariate file for training samples by subsetting eigenvec.
# PLINK expects covar file with FID IID then covariate columns.
# eigenvec has: FID IID PC1 PC2 ...

covar_file="$OUT_DIR/train_pca_covar.tsv"

# header
{
  printf 'FID\tIID'
  i=1
  while [[ $i -le $PCS ]]; do
    printf '\tPC%d' "$i"
    i=$(( i + 1 ))
  done
  printf '\n'
} > "$covar_file"

awk 'NR==FNR{key[$1"\t"$2]=1; next} {k=$1"\t"$2; if(k in key){
  printf "%s\t%s", $1, $2;
  for(i=1;i<='"$PCS"';i++) printf "\t%s", $(2+i);
  printf "\n";
}}' "$keep_train" "$PCA_EIGENVEC" >> "$covar_file"

# Build covar-name list
covar_names="PC1"
for ((i=2; i<=PCS; i++)); do
  covar_names+=",PC$i"
done

echo "[split_train_pca_gwas.sh] covar_names=${covar_names}"

# GWAS filters
filt_args=()
[[ -n "$MAF" && "$MAF" != "" && "$MAF" != "None" && "$MAF" != "NONE" ]] && filt_args+=(--maf "$MAF")
[[ -n "$GENO" && "$GENO" != "" && "$GENO" != "None" && "$GENO" != "NONE" ]] && filt_args+=(--geno "$GENO")
[[ -n "$MIND" && "$MIND" != "" && "$MIND" != "None" && "$MIND" != "NONE" ]] && filt_args+=(--mind "$MIND")
[[ -n "$HWE" && "$HWE" != "" && "$HWE" != "None" && "$HWE" != "NONE" ]] && filt_args+=(--hwe "$HWE")

# Run GWAS on training set using covariates from precomputed PCA
# Assumes phenotype is in .fam

gwas_prefix="$OUT_DIR/gwas"

$PLINK_CMD --bfile "$train_prefix" \
  --allow-no-sex \
  --covar "$covar_file" \
  --covar-name "$covar_names" \
  --logistic hide-covar \
  "${filt_args[@]}" \
  ${MEMORY:-} ${THREADS:-} \
  --out "$gwas_prefix"

mkdir -p "$(dirname "$DONE")"
echo "ok" > "$DONE"
