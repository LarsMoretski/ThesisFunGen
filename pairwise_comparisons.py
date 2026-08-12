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

## Import packages
from pathlib import Path
import numpy as np
import polars as pl
from tqdm import tqdm

## Hardcoded variables
parquet_dir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/kmerdatabase/10kreads_tryout")
output_dir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/kmerdatabase/outfiles")
output_dir.mkdir(parents=True,exist_ok=True)
min_count = 5 #Changes the amount of min reads per kmer to be allowed
## Define functions

def compare_kmers(df1, df2):
    '''
    Compares the presence of k-mers count tables (as dataframes)
    Returns the merged k-mer table for calculating similarities.
    '''

    merged = df1.join(df2, on="kmer", how="full", suffix="_2")
    merged = merged.fill_null(0)
    return merged

def calc_jaccard(merged):
    shared = ((pl.col("count") > 0) & (pl.col("count_2") > 0))
    union = merged.height
    intersection = merged.filter(shared).height
    if union == 0:
        return 0
    jaccard = intersection/union
    assert 0 <= jaccard <= 1, "Jaccard value must be between 0 and 1"
    return jaccard

def calc_weighted_jaccard(merged):
    numerator = (merged.select(pl.min_horizontal("count", "count_2").sum()).item())
    denominator = (merged.select(pl.max_horizontal("count", "count_2").sum()).item())
    if denominator == 0:
        return 0
    weighted_jaccard = numerator/denominator
    assert 0 <= weighted_jaccard <= 1, "Weighted Jaccard value must be between 0 and 1"
    return weighted_jaccard

def calc_dice_sorensen(merged):
    shared = ((pl.col("count") > 0) & (pl.col("count_2") > 0))
    intersection = merged.filter(shared).height
    sample1 = merged.filter(pl.col("count") > 0).height
    sample2 = merged.filter(pl.col("count_2") > 0).height
    denominator = sample1 + sample2
    if denominator ==0:
        return 0
    dice_coeff = (2*intersection)/denominator
    assert 0 <= dice_coeff <= 1, "Dice-Sørensen coefficient must be between 0 and 1"
    return dice_coeff

def write_to_matrix(matrix, sample_ids, outfile):
    df = pl.DataFrame(matrix, schema=sample_ids)
    df.insert_column(0, pl.Series("sample", sample_ids))
    df.write_csv(outfile, separator="\t")


## Main

def main():
    samples = {}

    for parquet_file in tqdm(list(parquet_dir.glob("*.parquet")), desc ="Loading k-mer database"):
        sample_id = parquet_file.stem
        samples[sample_id] = (pl.read_parquet(parquet_file).filter(pl.col("count") >= min_count))

    print(f"Loaded {len(samples)} samples")

    sample_ids = sorted(samples.keys())
    n_samples = len(sample_ids)

    jaccard_matrix = np.full((n_samples, n_samples), np.nan)
    weighted_jaccard_matrix = np.full((n_samples, n_samples), np.nan)
    dice_matrix = np.full((n_samples, n_samples), np.nan)
    np.fill_diagonal(jaccard_matrix, 1)
    np.fill_diagonal(weighted_jaccard_matrix, 1)
    np.fill_diagonal(dice_matrix, 1)

    total_comparisons = n_samples * (n_samples -1) // 2
    with tqdm(total=total_comparisons, desc= "Calculating k-mer similarities") as pbar:
        for i, sample1 in enumerate(sample_ids):
            for j, sample2 in enumerate(sample_ids[i+1:], start=i+1):
                df1 = samples[sample1]
                df2 = samples[sample2]
                merged = compare_kmers(df1, df2)

                jaccard_matrix[i,j]= calc_jaccard(merged)
                weighted_jaccard_matrix[i,j]= calc_weighted_jaccard(merged)
                dice_matrix[i,j]= calc_dice_sorensen(merged)
                pbar.update(1)

    write_to_matrix(jaccard_matrix, sample_ids, output_dir / "jaccard_matrix.tsv")
    write_to_matrix(weighted_jaccard_matrix, sample_ids, output_dir / "weighted_jaccard_matrix.tsv")
    write_to_matrix(dice_matrix, sample_ids, output_dir / "dice_matrix.tsv")

##
if __name__ == "__main__":
    main()

