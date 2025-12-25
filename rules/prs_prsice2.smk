import sys
from pathlib import Path

sys.path.append("utils")
from bioconfigme import (
    get_analysis_value,
    get_results_dir,
)


RESULTS_DIR = Path(get_results_dir())

DEFAULT_RESOURCES = get_analysis_value("default_resources", default={})
MEM_MB = int(DEFAULT_RESOURCES.get("mem_mb", 32000))
CORES = int(DEFAULT_RESOURCES.get("cores", 2))
TIME = str(DEFAULT_RESOURCES.get("time", "00:30:00"))

TARGETS_CFG = get_analysis_value("target_analysis", default={})
if not isinstance(TARGETS_CFG, dict) or not TARGETS_CFG:
    raise ValueError("analysis.yml must define target_analysis with at least one target")

PRSICE2_CFG = get_analysis_value("prsice2", default={})
if not isinstance(PRSICE2_CFG, dict):
    PRSICE2_CFG = {}

PRSICE_BIN = str(PRSICE2_CFG.get("prsice_bin", "bin/PRSice"))
PRSICE_R = str(PRSICE2_CFG.get("prsice_r", "bin/PRSice.R"))
RSCRIPT = str(PRSICE2_CFG.get("rscript", "Rscript"))
FORMATTER_R = str(PRSICE2_CFG.get("formatter_r", "scripts/format_gwas_for_prsice.R"))
PRSICE_MODULE = str(PRSICE2_CFG.get("module", ""))

STAT = str(PRSICE2_CFG.get("stat", "OR"))
SCORE = str(PRSICE2_CFG.get("score", "avg"))
BINARY_TARGET = str(PRSICE2_CFG.get("binary_target", "T"))

CLUMP_KB = str(PRSICE2_CFG.get("clump_kb", "250"))
CLUMP_R2 = str(PRSICE2_CFG.get("clump_r2", "0.1"))
CLUMP_P = str(PRSICE2_CFG.get("clump_p", "1"))
BAR_LEVELS = str(PRSICE2_CFG.get("bar_levels", ""))
PERM = str(PRSICE2_CFG.get("perm", ""))
THREAD = str(PRSICE2_CFG.get("thread", ""))
ADDITIONAL_PARAMS = str(PRSICE2_CFG.get("additional_params", ""))
QUANTILE = str(PRSICE2_CFG.get("quantile", ""))
QUANT_BREAK = str(PRSICE2_CFG.get("quant_break", ""))
QUANT_REF = str(PRSICE2_CFG.get("quant_ref", ""))


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


def _tools_for_target(target: str):
    cfg = TARGETS_CFG.get(target, {})
    if not isinstance(cfg, dict):
        return []
    tools = cfg.get("tools", [])
    if tools is None:
        return []
    if isinstance(tools, str):
        return [tools]
    if isinstance(tools, list):
        return [str(t) for t in tools]
    return []


def _ld_prefix_for_target(target: str) -> str:
    cfg = TARGETS_CFG.get(target, {})
    if not isinstance(cfg, dict):
        return ""
    sub = cfg.get("prsice2", {})
    if not isinstance(sub, dict):
        return ""
    return str(sub.get("ld_prefix", ""))


TARGETS = sorted([t for t in TARGETS_CFG.keys() if "PRSice2" in _tools_for_target(t)])

TOOL_DONE = [str(RESULTS_DIR / t / "prs" / "PRSice2" / "prs.done") for t in TARGETS]


rule all:
    input:
        TOOL_DONE


rule prsice2_bootstrap:
    input:
        gwas_assoc=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "gwas.assoc.logistic"
        ),
        test_bim=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "test.bim"
        ),
        test_bed=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "test.bed"
        ),
        test_fam=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "test.fam"
        ),
    output:
        done=str(
            RESULTS_DIR
            / "{target}"
            / "prs"
            / "PRSice2"
            / "bootstrap_{rep}"
            / "prsice2.done"
        )
    log:
        str(RESULTS_DIR / "{target}" / "log" / "prsice2_bootstrap_{rep}.log")
    params:
        out_dir=str(RESULTS_DIR / "{target}" / "prs" / "PRSice2" / "bootstrap_{rep}"),
        test_prefix=str(
            RESULTS_DIR
            / "{target}"
            / "split_train_pca_gwas"
            / "bootstrap_{rep}"
            / "test"
        ),
        prsice_bin=PRSICE_BIN,
        prsice_r=PRSICE_R,
        rscript=RSCRIPT,
        formatter_r=FORMATTER_R,
        prsice_module=PRSICE_MODULE,
        stat=STAT,
        score=SCORE,
        binary_target=BINARY_TARGET,
        clump_kb=CLUMP_KB,
        clump_r2=CLUMP_R2,
        clump_p=CLUMP_P,
        bar_levels=BAR_LEVELS,
        perm=PERM,
        thread=THREAD,
        quantile=QUANTILE,
        quant_break=QUANT_BREAK,
        quant_ref=QUANT_REF,
        additional_params=ADDITIONAL_PARAMS,
        ld_prefix=lambda wc: _ld_prefix_for_target(wc.target),
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

        echo "[prsice2_bootstrap] host=$(hostname)"
        echo "[prsice2_bootstrap] pwd=$(pwd)"
        echo "[prsice2_bootstrap] date=$(date)"

        out_dir="{params.out_dir}"
        mkdir -p "$out_dir"

        # Optional module load for R/PRSice environment
        if [[ -n "{params.prsice_module}" ]]; then
          if ! command -v module >/dev/null 2>&1; then
              [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
              [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
          fi
          module load {params.prsice_module}
        fi

        base_out="$out_dir/gwas.prsice.base.tsv"

                "{params.rscript}" "{params.formatter_r}" \
                    --in "{input.gwas_assoc}" \
                    --bim "{input.test_bim}" \
                    --out "$base_out"

        bash scripts/run_prsice2.sh \
          --prsice-bin "{params.prsice_bin}" \
          --prsice-r "{params.prsice_r}" \
          --rscript "{params.rscript}" \
          --base "$base_out" \
                    --target "{params.test_prefix}" \
          --out "$out_dir/prsice2" \
          --stat "{params.stat}" \
          --score "{params.score}" \
          --binary-target "{params.binary_target}" \
          --clump-kb "{params.clump_kb}" \
          --clump-r2 "{params.clump_r2}" \
          --clump-p "{params.clump_p}" \
          --bar-levels "{params.bar_levels}" \
          --perm "{params.perm}" \
          --thread "{params.thread}" \
                    --quantile "{params.quantile}" \
                    --quant-break "{params.quant_break}" \
                    --quant-ref "{params.quant_ref}" \
          --ld "{params.ld_prefix}" \
                    --extra "{params.additional_params}" \
          --done "{output.done}"
        """


rule prsice2_done:
    input:
        lambda wc: [
            str(
                RESULTS_DIR
                / wc.target
                / "prs"
                / "PRSice2"
                / f"bootstrap_{rep}"
                / "prsice2.done"
            )
            for rep in range(1, _n_bootstraps_for_target(wc.target) + 1)
        ]
    output:
        str(RESULTS_DIR / "{target}" / "prs" / "PRSice2" / "prs.done")
    log:
        str(RESULTS_DIR / "{target}" / "log" / "prsice2.done.log")
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        mkdir -p "$(dirname {log})"
        echo "[prsice2_done] date=$(date)" > {log}
        echo "ok" > {output}
        """
