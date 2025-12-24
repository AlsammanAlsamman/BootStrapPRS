#!/bin/bash
set -euo pipefail

CORES="${CORES:-2}"

# Runs the prune+pca+bootstrap targets from rules/prune_pca.smk.

snakemake \
  --snakefile rules/prune_pca.smk \
  --cores "$CORES" \
  --printshellcmds \
  --rerun-incomplete \
  --keep-going \
  all
