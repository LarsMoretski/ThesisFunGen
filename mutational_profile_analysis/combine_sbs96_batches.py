#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script for combining SBS96 mutational profiles of those run in batches
    isolates were run in batches for optimalization and checking inbetween runs

Usage: python3 combine_sbs96_batches.py
'''

## import packages
from pathlib import Path
import sys
import pandas as pd


## set path
ROOT = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
BATCH_ROOT = ROOT / "mutational_patterns" / "full_batches"
OUTPUT_DIR = ROOT / "mutational_patterns" / "combined"
OUTPUT_FILE = OUTPUT_DIR / "All_isolates.SBS96.tsv"


## def main
def main():
    matrices = []
    seen_samples = set()
    for batch_dir in sorted(BATCH_ROOT.glob("Batch[0-9]*")):
        project = batch_dir.name
        matrix_file = (batch_dir / "results" / "SBS" / f"{project}.SBS96.all")
        if not matrix_file.is_file():
            print(f"Missing matrix: {matrix_file}")
            sys.exit(1)

        df = pd.read_csv(matrix_file, sep="\t", index_col=0)
        if df.shape[0] != 96:
            print(f"Something has gone wrong, {project} has {df.shape[0]} contexts, not 96")
            sys.exit(1)

        duplicates = seen_samples.intersection(df.columns)
        if duplicates:
            print(f"Duplicate sample columns detected {sorted(duplicates)[:10]}")
            sys.exit(1)

        seen_samples.update(df.columns)
        matrices.append(df)
        print(f"Loaded {project}, {df.shape[1]} isolates")

    if not matrices:
        print("No matrices found.")
        sys.exit(1)

    reference_contexts = list(matrices[0].index) ## align rows by mutation context labels before combining
    for i, df in enumerate(matrices, start=1):
        if set(df.index) != set(reference_contexts):
            print(f"ERROR: Context labels differ in matrix {i}")
            sys.exit(1)

    combined = pd.concat(
        [df.loc[reference_contexts] for df in matrices], axis=1)

    if combined.shape != (96, 1260):
        print(f"The shape of the combined dimensions is not in the right format, {combined.shape}")

    if combined.columns.duplicated().any():
        print("Duplicate columns in combined matrix.")
        sys.exit(1)

    if (combined.select_dtypes(include="number") < 0).any().any():
        print("Negative mutation counts found")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    combined.to_csv(OUTPUT_FILE, sep="\t")

##
if __name__ == "__main__":
    main()
