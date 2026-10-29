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

## import packages
from pathlib import Path
import numpy as np
import pandas as pd

## set paths
BASE = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
OUT = BASE / "final_subdata/vcf_validation"
GENOTYPES_FILE = OUT / "genotypes.tsv"
SAMPLES_FILE = OUT / "snp_matrix_sample_order.txt"


samples = pd.read_csv(SAMPLES_FILE, header=None, dtype=str)[0].tolist()
n = len(samples)
differences = np.zeros((n, n), dtype=np.int64)
comparable = np.zeros((n, n), dtype=np.int64)
n_records = 0
n_snp_records = 0

with open(GENOTYPES_FILE) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        chrom, pos = fields[:2]
        gt = fields[2:]

        if len(gt) != n:
            raise ValueError(f"{chrom}:{pos}: expected {n} genotypes found {len(gt)} ")

        n_records += 1
        values = np.full(n, -1, dtype=np.int8)
        for i, g in enumerate(gt):
            if g in ("0", "1"):
                values[i] = int(g)
            elif g in ("0/0", "1/1", "0|0", "1|1"):
                values[i] = int(g[0])
            elif g in ("0/1", "1/0", "0|1", "1|0"):
                raise ValueError(f"Unexpected heterozygous genotype {g} at {chrom}:{pos}, sample {samples[i]} ")

        called_idx = np.flatnonzero(values >= 0)
        if len(called_idx) < 2:
            continue

        called_values = values[called_idx]
        comparable[np.ix_(called_idx, called_idx)] += 1
        differences[np.ix_(called_idx, called_idx)] += (called_values[:, None] != called_values[None, :])
        n_snp_records += 1

dissimilarity = np.full((n, n), np.nan, dtype=float)
np.divide(differences, comparable, out=dissimilarity, where=comparable > 0)
np.fill_diagonal(dissimilarity, 0.0)

for name, matrix in [("snp_differences", differences), ("comparable_sites", comparable), ("snp_dissimilarity", dissimilarity),]:
    df = pd.DataFrame(matrix, index=samples, columns=samples)
    df.index.name = "sample_id"
    df.to_csv(OUT / f"{name}.tsv", sep="\t", na_rep="NA")

pd.Series(samples, name="vcf_id").to_csv(OUT / "snp_matrix_sample_order.txt", index=False, header=False)
upper = np.triu_indices(n, 1)

print("Median comparable sites per pair:", int(np.median(comparable[upper])))
print("Median SNP differences per pair:", float(np.median(differences[upper])))
print("Median SNP dissimilarity:", float(np.nanmedian(dissimilarity[upper])))

