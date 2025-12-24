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

# Global bootstrap defaults (optional)
TRAIN_FRAC = float(get_analysis_value("bootstrap.train_frac", default=0.8))
SEED_BASE = int(get_analysis_value("bootstrap.seed_base", default=12345))

# Phenotype config (needed for GWAS/PRS steps)
PHENO_FILE = str(get_analysis_value("phenotype.file", default="inputs/pheno.tsv"))
PHENO_COL = str(get_analysis_value("phenotype.col", default="PHENO"))

# LD pruning parameters (optional)
PRUNE_WINDOW = int(get_analysis_value("prune.window_kb", default=200))
PRUNE_STEP = int(get_analysis_value("prune.step", default=50))
PRUNE_R2 = float(get_analysis_value("prune.r2", default=0.2))


def _bfile_inputs(prefix: str):
    return [f"{prefix}.bed", f"{prefix}.bim", f"{prefix}.fam"]


def _target_results_dir(target: str) -> Path:
    return RESULTS_DIR / target


def _target_log_dir(target: str) -> Path:
    return _target_results_dir(target) / "log"


def _target_qc_dir(target: str) -> Path:
    return _target_results_dir(target) / "qc"


def _target_bootstrap_dir(target: str) -> Path:
    return _target_results_dir(target) / "bootstrap"


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


PLINK_MODULE = get_software_module("plink2")
PLINK_CMD = get_software_command("plink2", default="plink")
PLINK_PARAMS = get_software_params("plink2")
PLINK_MEMORY = str(PLINK_PARAMS.get("memory", ""))
PLINK_THREADS = str(PLINK_PARAMS.get("threads", ""))


ALL_DONE = []
for _t in TARGETS:
    for _rep in range(1, _n_bootstraps_for_target(_t) + 1):
        ALL_DONE.append(str(_target_bootstrap_dir(_t) / f"rep={_rep}" / "bootstrap.done"))


rule all:
    input:
        ALL_DONE


rule prune_pca_done:
    input:
        prune_in=str(RESULTS_DIR / "{target}" / "prune_pca" / "prune.prune.in"),
        eigenvec=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenvec"),
        eigenval=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenval")
    output:
        done=str(RESULTS_DIR / "{target}" / "prune_pca" / "prune_pca.done")
    log:
        str(RESULTS_DIR / "{target}" / "prune_pca" / "prune_pca.log")
    resources:
        mem_mb=MEM_MB,
        time=TIME,
        cores=CORES
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done})" "$(dirname {log})"
        echo "ok" > {output.done}
        """  


rule prune_ld:
    input:
        bed=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bed",
        bim=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bim",
        fam=lambda wc: f"{_plink_prefix_for_target(wc.target)}.fam"
    output:
        prune_in=str(RESULTS_DIR / "{target}" / "prune_pca" / "prune.prune.in"),
        prune_out=str(RESULTS_DIR / "{target}" / "prune_pca" / "prune.prune.out"),
        pruned_bed=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.bed"),
        pruned_bim=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.bim"),
        pruned_fam=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.fam")
    log:
        str(RESULTS_DIR / "{target}" / "prune_pca" / "prune_ld.log")
    params:
        step_dir=str(RESULTS_DIR / "{target}" / "prune_pca")
    resources:
        mem_mb=MEM_MB,
        time=TIME,
        cores=CORES
    threads: CORES
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        mkdir -p {params.step_dir}

        bfile="{input.bed}"
        bfile="${{bfile%.bed}}"

        sorted_prefix="{params.step_dir}/sorted"

        if ! command -v module >/dev/null 2>&1; then
            [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
            [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
        fi
        module load {PLINK_MODULE}

        # PLINK LD pruning requires a sorted .bim. Create a sorted bed set first.
        {PLINK_CMD} --bfile "$bfile" \
            --make-bed \
            {PLINK_MEMORY} {PLINK_THREADS} \
            --out "$sorted_prefix"

        {PLINK_CMD} --bfile "$sorted_prefix" \
            --indep-pairwise {PRUNE_WINDOW} {PRUNE_STEP} {PRUNE_R2} \
            {PLINK_MEMORY} {PLINK_THREADS} \
            --out {params.step_dir}/prune

        {PLINK_CMD} --bfile "$sorted_prefix" \
            --extract {params.step_dir}/prune.prune.in \
            --make-bed \
            {PLINK_MEMORY} {PLINK_THREADS} \
            --out {params.step_dir}/pruned
        """


rule pca_once:
    input:
        bed=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.bed"),
        bim=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.bim"),
        fam=str(RESULTS_DIR / "{target}" / "prune_pca" / "pruned.fam")
    output:
        eigenvec=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenvec"),
        eigenval=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenval")
    log:
        str(RESULTS_DIR / "{target}" / "prune_pca" / "pca_once.log")
    params:
        pcs=lambda wc: _pcs_for_target(wc.target),
        step_dir=str(RESULTS_DIR / "{target}" / "prune_pca")
    resources:
        mem_mb=MEM_MB,
        time=TIME,
        cores=CORES
    threads: CORES
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        mkdir -p {params.step_dir}

        if ! command -v module >/dev/null 2>&1; then
            [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
            [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
        fi
        module load {PLINK_MODULE}

        {PLINK_CMD} --bfile {params.step_dir}/pruned \
            --pca {params.pcs} \
            {PLINK_MEMORY} {PLINK_THREADS} \
            --out {params.step_dir}/pca
        """


rule bootstrap_gwas_prs_rep:
    input:
        pca_eigenvec=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenvec"),
        pca_eigenval=str(RESULTS_DIR / "{target}" / "prune_pca" / "pca.eigenval"),
        bed=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bed",
        bim=lambda wc: f"{_plink_prefix_for_target(wc.target)}.bim",
        fam=lambda wc: f"{_plink_prefix_for_target(wc.target)}.fam"
    output:
        done=str(RESULTS_DIR / "{target}" / "bootstrap" / "rep={rep}" / "bootstrap.done")
    log:
        str(RESULTS_DIR / "{target}" / "log" / "bootstrap_rep={rep}.log")
    params:
        train_frac=TRAIN_FRAC,
        pcs=lambda wc: _pcs_for_target(wc.target),
        log_dir=str(RESULTS_DIR / "{target}" / "log")
    resources:
        mem_mb=MEM_MB,
        time=TIME,
        cores=CORES
    threads: CORES
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        mkdir -p {params.log_dir}

        bfile="{input.bed}"
        bfile="${{bfile%.bed}}"

        out_dir="$(dirname {output.done})"

        seed=$(( {SEED_BASE} + {wildcards.rep} ))

        if ! command -v module >/dev/null 2>&1; then
            [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
            [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
        fi
        module load {PLINK_MODULE}

        python scripts/example_helper.py \
            --bfile "$bfile" \
            --pca-eigenvec {input.pca_eigenvec} \
            --pheno-file {PHENO_FILE} \
            --pheno-col {PHENO_COL} \
            --train-frac {params.train_frac} \
            --seed "$seed" \
            --max-pcs {params.pcs} \
            --plink-cmd {PLINK_CMD} \
            --out-dir "$out_dir" \
            --done {output.done}
        """
