#!/bin/bash
# Convenience runner for PRSice2 step (produces the tool-level done marker).
set -euo pipefail

./submit.sh --snakefile rules/prs_prsice2.smk results/LAMR_PRS/prs/PRSice2/prs.done "$@"
