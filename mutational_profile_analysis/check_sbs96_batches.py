#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 30/07/2026

Description:
    Python script for computing the pairwise distance between two k-mer indexes using:
        - Jaccard distance
        - Weighted Jaccard similarity
        - Dice-Sørensen coefficient

Usage: python3 jellyfish_to_parquet.py --indir /path/to/files.jf --outdir /path/to/outdir
'''

## import statements
from pathlib import Path
import sys
import pandas as pd

## Set paths
ROOT = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
BATCH_ROOT = ROOT / "mutational_patterns" / "full_batches"
EXPECTED_CONTEXTS = 96


## main
def main():
    total_samples = 0
    total_expected = 0
    problems = []
    seen_samples = {}
    batch_dirs = sorted(BATCH_ROOT.glob("Batch[0-9]*"))

    if not batch_dirs:
        print(f"ERROR: No batch directories found in {BATCH_ROOT}")
        sys.exit(1)

    print("SigProfiler SBS96 batch validation")
    print("=" * 60)

    for batch_dir in batch_dirs:
        project = batch_dir.name
        matrix_file = (batch_dir / "results" / "SBS" / f"{project}.SBS96.all")
        manifest_file = batch_dir / "vcfs.txt"
        complete_file = batch_dir / "COMPLETE"

        if not complete_file.is_file():
            problems.append(f"{project}: missing COMPLETE marker")

        if not matrix_file.is_file() or matrix_file.stat().st_size == 0:
            print(f"FAIL  {project}: missing or empty matrix")
            continue

        if not manifest_file.is_file():
            problems.append(f"{project}: missing vcfs.txt manifest")
            expected_samples = None
        else:
            with manifest_file.open() as handle:
                expected_samples = sum(
                    1 for line in handle if line.strip()
                )
            total_expected += expected_samples

        try:
            df = pd.read_csv(matrix_file, sep="\t", index_col=0)
        except Exception as exc:
            problems.append(f"{project}: matrix read error: {exc}")
            print(f"FAIL  {project}: could not read matrix: {exc}")
            continue

        n_contexts, n_samples = df.shape
        total_samples += n_samples

        status = "OK"
        if n_contexts != EXPECTED_CONTEXTS:
            problems.append(f"{project}: expected {EXPECTED_CONTEXTS} contexts found i {n_contexts} ")
            status = "FAIL"

        if expected_samples is not None and n_samples != expected_samples:
            problems.append(f"{project}: expected {expected_samples} sample columns found in {n_samples} ")
            status = "FAIL"

        if df.columns.duplicated().any():
            problems.append(f"{project}: duplicate sample columns" )
            status = "FAIL"

        for sample in df.columns:
            if sample in seen_samples:
                problems.append(f"Duplicate sample across batches {sample} {seen_samples[sample]} and {project} ")
                status = "FAIL"
            else:
                seen_samples[sample] = project

    print("=" * 60)
    print(f"Batches found: {len(batch_dirs)}")
    print(f"Expected sample columns: {total_expected}")
    print(f"Observed sample columns: {total_samples}")
    print(f"Unique sample names: {len(seen_samples)}")

    if problems:
        print(f"\nValidation found {len(problems)} issue(s):")
        for problem in problems:
            print(f" - {problem}")
        sys.exit(1)

    print("All batch matrices passed the checks.")

## main
if __name__ == "__main__":
    main()
