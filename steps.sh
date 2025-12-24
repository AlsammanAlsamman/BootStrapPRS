#!/bin/bash
# Step 1: Prune + PCA (target: LAMR_PRS)
./submit.sh --snakefile rules/prune_pca.smk results/LAMR_PRS/prune_pca/prune_pca.done

# Step 2: Split train/test + GWAS on training (all bootstraps)
./submit.sh --snakefile rules/split_train_pca_gwas.smk results/LAMR_PRS/split_train_pca_gwas/split_train_pca_gwas.done --jobs 20

# Step 3: PRS (PRSice2) on test set (all bootstraps)
./submit.sh --snakefile rules/prs_prsice2.smk results/LAMR_PRS/prs/PRSice2/prs.done --jobs 20

