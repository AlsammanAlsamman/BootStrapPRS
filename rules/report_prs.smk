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

RSCRIPT = str(PRSICE2_CFG.get("rscript", "Rscript"))
PRSICE_MODULE = str(PRSICE2_CFG.get("module", ""))


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


TARGETS_PRSICE2 = sorted([t for t in TARGETS_CFG.keys() if "PRSice2" in _tools_for_target(t)])

REPORT_DONE = [
    str(RESULTS_DIR / t / "prs" / "PRSice2" / "report" / "report.done")
    for t in TARGETS_PRSICE2
]


rule all:
    input:
        REPORT_DONE


rule report_prsice2:
    input:
        split_done=str(RESULTS_DIR / "{target}" / "split_train_pca_gwas" / "split_train_pca_gwas.done"),
    output:
        done=str(RESULTS_DIR / "{target}" / "prs" / "PRSice2" / "report" / "report.done"),
    log:
        str(RESULTS_DIR / "{target}" / "log" / "report_prsice2.log"),
    params:
        prsice_module=PRSICE_MODULE,
        rscript=RSCRIPT,
        results_dir=str(RESULTS_DIR),
        prsice2_dir=str(RESULTS_DIR / "{target}" / "prs" / "PRSice2"),
        split_dir=str(RESULTS_DIR / "{target}" / "split_train_pca_gwas"),
        out_dir=str(RESULTS_DIR / "{target}" / "prs" / "PRSice2" / "report"),
        bootstraps=lambda wc: _n_bootstraps_for_target(wc.target),
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

        echo "[report_prsice2] host=$(hostname)"
        echo "[report_prsice2] pwd=$(pwd)"
        echo "[report_prsice2] date=$(date)"

        mkdir -p "{params.out_dir}"

        # Optional module load for R environment
        if [[ -n "{params.prsice_module}" ]]; then
          if ! command -v module >/dev/null 2>&1; then
              [ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh
              [ -f /usr/share/Modules/init/bash ] && . /usr/share/Modules/init/bash
          fi
          module load {params.prsice_module}
        fi

        "{params.rscript}" scripts/report_prsice2_bootstraps.R \
                    --results-dir "{params.results_dir}" \
          --target "{wildcards.target}" \
          --bootstraps "{params.bootstraps}" \
          --prsice2-dir "{params.prsice2_dir}" \
          --split-dir "{params.split_dir}" \
          --out-dir "{params.out_dir}"

        mkdir -p "$(dirname {output.done})"
        echo "ok" > "{output.done}"
        """
