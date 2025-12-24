#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF' >&2
Format PLINK .assoc.logistic to PRSice2 base format with required columns.

Usage:
  bash format_plink_gwas_for_prsice.sh --assoc GWAS.assoc.logistic --bim target.bim --out base.tsv

Output columns (tab-delimited):
  SNP  CHR  BP  A1  A2  OR  P

Notes:
- Keeps only TEST==ADD rows.
- A2 is pulled from the BIM file.
EOF
  exit 2
}

assoc=""
bim=""
out=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --assoc) assoc="$2"; shift 2;;
    --bim) bim="$2"; shift 2;;
    --out) out="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1" >&2; usage;;
  esac
done

[[ -n "$assoc" && -n "$bim" && -n "$out" ]] || usage
[[ -f "$assoc" ]] || { echo "Missing assoc: $assoc" >&2; exit 1; }
[[ -f "$bim" ]] || { echo "Missing bim: $bim" >&2; exit 1; }

mkdir -p "$(dirname "$out")"

awk -v bim="$bim" '
BEGIN{
  FS=OFS="\t";
  # Load SNP->A2 from bim (col2=snp, col6=a2)
  while ((getline line < bim) > 0) {
    n=split(line, a, /[ \t]+/);
    if (n>=6) a2[a[2]]=a[6];
  }
  close(bim);
}
NR==1{
  # detect header indices in assoc
  FS=OFS="\t";
  for(i=1;i<=NF;i++) {
    h=$i;
    idx[h]=i;
  }
  # required columns
  if (!("SNP" in idx) || !("CHR" in idx) || !("BP" in idx) || !("A1" in idx) || !("TEST" in idx) || !("OR" in idx) || !("P" in idx)) {
    print "ERROR: assoc header missing required columns" > "/dev/stderr";
    print "Have: " $0 > "/dev/stderr";
    exit 2;
  }
  print "SNP","CHR","BP","A1","A2","OR","P";
  next;
}
{
  test=$idx["TEST"];
  if (test != "ADD") next;
  snp=$idx["SNP"]; chr=$idx["CHR"]; bp=$idx["BP"]; a1=$idx["A1"]; orv=$idx["OR"]; pv=$idx["P"];
  if (snp=="" || pv=="NA") next;
  a2v = (snp in a2 ? a2[snp] : "NA");
  print snp, chr, bp, a1, a2v, orv, pv;
}
' "$assoc" > "$out"
