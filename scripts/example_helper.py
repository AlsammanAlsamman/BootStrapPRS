#!/usr/bin/env python

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import random
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


def _run(cmd: Sequence[str], *, cwd: Optional[str] = None) -> None:
    subprocess.run(list(cmd), cwd=cwd, check=True)


def _read_fam_iids(fam_path: Path) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    with fam_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            pairs.append((parts[0], parts[1]))
    if not pairs:
        raise RuntimeError(f"No samples found in {fam_path}")
    return pairs


def _write_keep(path: Path, pairs: Iterable[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for fid, iid in pairs:
            f.write(f"{fid}\t{iid}\n")


def _read_eigenvec(path: Path) -> Tuple[List[str], Dict[Tuple[str, str], List[str]]]:
    """Returns (pc_names, map[(FID,IID)] -> [pc1..pcK])."""
    sample_to_pcs: Dict[Tuple[str, str], List[str]] = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            fid, iid = parts[0], parts[1]
            pcs = parts[2:]
            sample_to_pcs[(fid, iid)] = pcs

    # Determine number of PCs from first entry
    if not sample_to_pcs:
        raise RuntimeError(f"No PCA rows found in {path}")
    first_pcs = next(iter(sample_to_pcs.values()))
    pc_names = [f"PC{i+1}" for i in range(len(first_pcs))]
    return pc_names, sample_to_pcs


def _write_covar(
    out_path: Path,
    keep_set: Set[Tuple[str, str]],
    pc_names: List[str],
    sample_to_pcs: Dict[Tuple[str, str], List[str]],
    max_pcs: int = 10,
) -> List[str]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    used_pc_names = pc_names[: min(max_pcs, len(pc_names))]

    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["FID", "IID", *used_pc_names])
        for fid, iid in sorted(keep_set):
            pcs = sample_to_pcs.get((fid, iid))
            if not pcs:
                continue
            writer.writerow([fid, iid, *pcs[: len(used_pc_names)]])

    return used_pc_names


def _find_glm_file(prefix: Path) -> Path:
    matches = glob.glob(str(prefix) + ".*.glm.*")
    if not matches:
        raise RuntimeError(f"No PLINK2 --glm output found for prefix: {prefix}")
    matches_sorted = sorted(matches)
    return Path(matches_sorted[0])


def _find_plink_assoc_logistic(prefix: Path) -> Path:
    candidates = [
        Path(str(prefix) + ".assoc.logistic"),
        Path(str(prefix) + ".assoc.logistic.gz"),
    ]
    for p in candidates:
        if p.exists():
            return p

    matches = glob.glob(str(prefix) + ".*.assoc.logistic*")
    if matches:
        return Path(sorted(matches)[0])
    raise RuntimeError(f"No PLINK --logistic output found for prefix: {prefix}")


def _write_score_from_glm(glm_path: Path, score_path: Path) -> None:
    """Create a 3-column score file: ID A1 BETA.

    Handles either BETA or OR columns (beta=log(OR)). Filters to TEST==ADD when available.
    """

    score_path.parent.mkdir(parents=True, exist_ok=True)

    with glm_path.open("r", encoding="utf-8") as f_in, score_path.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        reader = csv.reader(f_in, delimiter="\t")
        writer = csv.writer(f_out, delimiter="\t")

        header: Optional[List[str]] = None
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                # Some PLINK2 outputs prefix header with '#'
                row[0] = row[0].lstrip("#")
            header = row
            break

        if header is None:
            raise RuntimeError(f"Empty glm file: {glm_path}")

        col = {name: idx for idx, name in enumerate(header)}
        required_any = ("BETA" in col) or ("OR" in col)
        if "ID" not in col or "A1" not in col or not required_any:
            raise RuntimeError(
                f"Unexpected glm columns in {glm_path}. Need ID, A1 and BETA or OR. Got: {header}"
            )

        test_idx = col.get("TEST")
        beta_idx = col.get("BETA")
        or_idx = col.get("OR")

        writer.writerow(["ID", "A1", "BETA"])

        for row in reader:
            if not row:
                continue
            if test_idx is not None and len(row) > test_idx:
                if row[test_idx] and row[test_idx] != "ADD":
                    continue

            variant_id = row[col["ID"]]
            a1 = row[col["A1"]]
            beta: float
            if beta_idx is not None and len(row) > beta_idx and row[beta_idx] not in ("", "NA"):
                beta = float(row[beta_idx])
            elif or_idx is not None and len(row) > or_idx and row[or_idx] not in ("", "NA"):
                beta = math.log(float(row[or_idx]))
            else:
                continue

            if not variant_id or not a1:
                continue
            writer.writerow([variant_id, a1, f"{beta}"])


def _write_score_from_assoc_logistic(assoc_path: Path, score_path: Path) -> None:
    """Create a 3-column score file from PLINK 1.9 .assoc.logistic: SNP A1 BETA.

    Uses log(OR) as BETA. Filters to TEST==ADD when present.
    """

    score_path.parent.mkdir(parents=True, exist_ok=True)

    with assoc_path.open("r", encoding="utf-8") as f_in, score_path.open(
        "w", encoding="utf-8", newline=""
    ) as f_out:
        header: Optional[List[str]] = None
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            header = line.split()
            break
        if header is None:
            raise RuntimeError(f"Empty assoc file: {assoc_path}")

        col = {name: idx for idx, name in enumerate(header)}
        if "SNP" not in col or "A1" not in col or "OR" not in col:
            raise RuntimeError(
                f"Unexpected assoc columns in {assoc_path}. Need SNP, A1, OR. Got: {header}"
            )

        test_idx = col.get("TEST")

        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["ID", "A1", "BETA"])

        for line in f_in:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if test_idx is not None and len(parts) > test_idx:
                if parts[test_idx] and parts[test_idx] != "ADD":
                    continue
            snp = parts[col["SNP"]]
            a1 = parts[col["A1"]]
            or_str = parts[col["OR"]]
            if or_str in ("NA", "nan", "NaN", ""):
                continue
            beta = math.log(float(or_str))
            if not snp or not a1:
                continue
            writer.writerow([snp, a1, f"{beta}"])


def main() -> int:
    ap = argparse.ArgumentParser(description="Bootstrap: split->GWAS->PRS for one replicate")
    ap.add_argument("--bfile", required=True, help="PLINK bed/bim/fam prefix")
    ap.add_argument("--pca-eigenvec", required=True, help="PCA eigenvec file from pruned data")
    ap.add_argument("--pheno-file", required=True, help="Phenotype file")
    ap.add_argument("--pheno-col", required=True, help="Phenotype column name")
    ap.add_argument("--train-frac", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out-dir", required=True, help="Output directory for this replicate")
    ap.add_argument("--done", required=True, help="Path to .done marker to create")
    ap.add_argument("--max-pcs", type=int, default=10, help="Number of PCs to include as covariates")
    ap.add_argument("--plink-cmd", default="plink", help="plink executable (plink or plink2)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fam_path = Path(f"{args.bfile}.fam")
    samples = _read_fam_iids(fam_path)

    rng = random.Random(args.seed)
    rng.shuffle(samples)

    n_train = max(1, int(round(len(samples) * args.train_frac)))
    n_train = min(n_train, len(samples) - 1) if len(samples) > 1 else 1

    train_pairs = samples[:n_train]
    test_pairs = samples[n_train:]

    keep_train = out_dir / "split" / "train.keep"
    keep_test = out_dir / "split" / "test.keep"
    _write_keep(keep_train, train_pairs)
    _write_keep(keep_test, test_pairs)

    pc_names, sample_to_pcs = _read_eigenvec(Path(args.pca_eigenvec))
    train_set = set(train_pairs)
    test_set = set(test_pairs)

    covar_train = out_dir / "covar" / "train_pca.tsv"
    covar_test = out_dir / "covar" / "test_pca.tsv"

    covar_cols = _write_covar(
        covar_train, train_set, pc_names, sample_to_pcs, max_pcs=args.max_pcs
    )
    _write_covar(covar_test, test_set, pc_names, sample_to_pcs, max_pcs=args.max_pcs)

    covar_name_arg = ",".join(covar_cols)

    gwas_prefix = out_dir / "gwas" / "train"
    gwas_prefix.parent.mkdir(parents=True, exist_ok=True)

    plink_cmd_lower = Path(args.plink_cmd).name.lower()
    is_plink2 = "plink2" in plink_cmd_lower

    if is_plink2:
        _run(
            [
                args.plink_cmd,
                "--bfile",
                args.bfile,
                "--keep",
                str(keep_train),
                "--pheno",
                args.pheno_file,
                "--pheno-name",
                args.pheno_col,
                "--covar",
                str(covar_train),
                "--covar-name",
                covar_name_arg,
                "--glm",
                "hide-covar",
                "cols=+a1freq",
                "--out",
                str(gwas_prefix),
            ]
        )
        glm_path = _find_glm_file(gwas_prefix)
        score_file = out_dir / "prs" / "scorefile.tsv"
        _write_score_from_glm(glm_path, score_file)
    else:
        _run(
            [
                args.plink_cmd,
                "--bfile",
                args.bfile,
                "--keep",
                str(keep_train),
                "--pheno",
                args.pheno_file,
                "--pheno-name",
                args.pheno_col,
                "--covar",
                str(covar_train),
                "--covar-name",
                covar_name_arg,
                "--logistic",
                "hide-covar",
                "--out",
                str(gwas_prefix),
            ]
        )
        assoc_path = _find_plink_assoc_logistic(gwas_prefix)
        score_file = out_dir / "prs" / "scorefile.tsv"
        _write_score_from_assoc_logistic(assoc_path, score_file)

    prs_prefix = out_dir / "prs" / "test"
    prs_prefix.parent.mkdir(parents=True, exist_ok=True)

    if is_plink2:
        _run(
            [
                args.plink_cmd,
                "--bfile",
                args.bfile,
                "--keep",
                str(keep_test),
                "--score",
                str(score_file),
                "1",
                "2",
                "3",
                "header-read",
                "cols=+scoresums",
                "--out",
                str(prs_prefix),
            ]
        )
    else:
        _run(
            [
                args.plink_cmd,
                "--bfile",
                args.bfile,
                "--keep",
                str(keep_test),
                "--score",
                str(score_file),
                "1",
                "2",
                "3",
                "header-read",
                "sum",
                "--out",
                str(prs_prefix),
            ]
        )

    done_path = Path(args.done)
    done_path.parent.mkdir(parents=True, exist_ok=True)
    done_path.write_text(
        f"ok\nseed={args.seed}\ntrain_n={len(train_pairs)}\ntest_n={len(test_pairs)}\n",
        encoding="utf-8",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
