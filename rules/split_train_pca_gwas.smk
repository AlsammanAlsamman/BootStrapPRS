import sys
from pathlib import Path

sys.path.append("utils")
from bioconfigme import (
    get_analysis_value,
    get_results_dir,
    get_software_command,
    get_software_module,
    get_software_params,
)


RESULTS_DIR = Path(get_results_dir())

DEFAULT_RESOURCES = get_analysis_value("default_resources", default={})
MEM_MB = int(DEFAULT_RESOURCES.get("mem_mb", 32000))
CORES = int(DEFAULT_RESOURCES.get("cores", 2))
TIME = str(DEFAULT_RESOURCES.get("time", "00:30:00"))

TARGETS_CFG = get_analysis_value("target_analysis", default={})
if not isinstance(TARGETS_CFG, dict) or not TARGETS_CFG:
    raise ValueError("analysis.yml must define target_analysis with at least one target")

TARGETS = sorted(TARGETS_CFG.keys())

BOOTSTRAP_CFG = get_analysis_value("bootstrap", default={})
if not isinstance(BOOTSTRAP_CFG, dict):
    BOOTSTRAP_CFG = {}

TRAIN_FRAC = float(BOOTSTRAP_CFG.get("train_frac", 0.8))
SEED_BASE = int(BOOTSTRAP_CFG.get("seed_base", 12345))

GWAS_FILTERS = get_analysis_value("gwas_filters", default={})
if not isinstance(GWAS_FILTERS, dict):
    GWAS_FILTERS = {}

MAF = str(GWAS_FILTERS.get("maf", ""))
GENO = str(GWAS_FILTERS.get("geno", ""))
MIND = str(GWAS_FILTERS.get("mind", ""))
HWE = str(GWAS_FILTERS.get("hwe", ""))

PLINK_MODULE = get_software_module("plink2")
PLINK_CMD = get_software_command("plink2", default="plink")
PLINK_PARAMS = get_software_params("plink2")
PLINK_MEMORY = str(PLINK_PARAMS.get("memory", ""))
PLINK_THREADS = str(PLINK_PARAMS.get("threads", ""))


def _plink_prefix_for_target(target: str) -> str:
    cfg = TARGETS_CFG.get(target, {})
    if not isinstance(cfg, dict) or "plink_prefix" not in cfg:
        raise KeyError(f"target_analysis.{target}.plink_prefix is required")
    return str(cfg["plink_prefix"])


def _n_bootstraps_for_target(target: str) -> int:
    cfg = TARGETS_CFG.get(target, {})
    if not isinstance(cfg, dict):
        return 1
    return int(cfg.get("bootstraps", 1))


def _pcs_for_target(target: str) -> int:
    cfg = TARGETS_CFG.get(target, {})
    if not isinstance(cfg, dict):
        return 10
    return int(cfg.get("PCs", 10))


TARGET_DONE = [
    str(RESULTS_DIR / _t / "split_train_pca_gwas" / "split_train_pca_gwas.done")
    for _t in TARGETS
]


rule all:
    input:
        TARGET_DONE


rule split_train_pca_gwas_done:
    input:
        lambda wc: [
            str(
                RESULTS_DIR
                / wc.target
                / "split_train_pca_gwas"
                / f"bootstrap_{rep}"
                / "split_train_pca_gwas.done"
            )
            for rep in range(1, _n_bootstraps_for_target(wc.target) + 1)
        ]
    output:
        str(RESULTS_DIR / "{target}" / "split_train_pca_gwas" / "split_train_pca_gwas.done")
    log:
        str(RESULTS_DIR / "{target}" / "split_train_pca_gwas" / "split_train_pca_gwas.done.log")
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        echo "[split_train_pca_gwas_done] date=$(date)" > {log}
        echo "ok" > {output}
        """


rule split_train_pca_gwas:
    input:
        pca_eigenvec=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenvec"),
        bed=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bed",
        bim=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bim",
        fam=lambda wc: f"{_plink_prefix_for_target(wc.target)}.fam",
    output:
        done=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "split_train_pca_gwas.done"
        )
    log:
        str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "split_train_pca_gwas.log"
        )
    params:
        pcs=lambda wc: _pcs_for_target(wc.target),
        bfile_prefix=lambda wc: _plink_prefix_for_target(wc.target),
        train_frac=TRAIN_FRAC,
        seed_base=SEED_BASE,
        maf=MAF,
        geno=GENO,
        mind=MIND,
        hwe=HWE,
        memory=PLINK_MEMORY,
        threads=PLINK_THREADS,
    resources:
        mem_mb=MEM_MB,
        time=TIME,
        cores=CORES
    threads: CORES
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        exec > {log} 2>&1
        set -euo pipefail

        echo "[split_train_pca_gwas] host=$(hostname)"
        echo "[split_train_pca_gwas] pwd=$(pwd)"
        echo "[split_train_pca_gwas] date=$(date)"

        out_dir="$(dirname {output.done})"
        mkdir -p "$out_dir"

        if ! command -v module >/dev/null 2>&1; then
            [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
            [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
        fi
        module load {PLINK_MODULE}

        bash scripts/split_train_pca_gwas.sh \
            --plink-cmd "{PLINK_CMD}" \
            --bfile "{params.bfile_prefix}" \
            --pca-eigenvec "{input.pca_eigenvec}" \
            --out-dir "$out_dir" \
            --rep "{wildcards.rep}" \
            --seed-base "{params.seed_base}" \
            --train-frac "{params.train_frac}" \
            --pcs "{params.pcs}" \
            --maf "{params.maf}" \
            --geno "{params.geno}" \
            --mind "{params.mind}" \
            --hwe "{params.hwe}" \
            --memory "{params.memory}" \
            --threads "{params.threads}" \
            --done "{output.done}"
        """
