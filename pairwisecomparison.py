#!/usr/bin/env python3

"""
Author: Lars Moret, MSc student @ FunGen

Description:
    Calculate pairwise k-mer similarities.

    Metrics:
        - Jaccard
        - Weighted Jaccard
        - Dice-Sørensen

    Only the upper triangle is calculated/written.
    The lower triangle is completed later by run_pipeline.py.

Resume behavior:
    - If all three matrices exist, the comparison is skipped.
    - If one or more are missing, the complete comparison is recalculated.

Usage:
    python3 pairwisecomparison.py \
        --indir /path/to/parquet \
        --outdir /path/to/matrices \
        --prefix 500k_27 \
        --threads 16 \
        --min-count 5
"""

## import packges
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import polars as pl
from tqdm import tqdm

SAMPLES = {}

## set output paths
def matrix_paths(outdir, prefix):
    return [
        outdir / f"{prefix}_jaccard.tsv",
        outdir / f"{prefix}_weighted_jaccard.tsv",
        outdir / f"{prefix}_dice.tsv"]


## Define functions
def compare_kmers(df1, df2):
    merged = df1.join(df2, on="kmer", how="full", suffix="_2")

    return merged.fill_null(0)



def calc_jaccard(merged):
    shared = ((pl.col("count") > 0) & (pl.col("count_2") > 0))
    union = merged.height
    intersection = (merged.filter(shared).height)
    if union == 0:
        return 0.0

    return intersection / union


def calc_weighted_jaccard(merged):
    numerator = (merged.select(pl.min_horizontal("count", "count_2").sum()).item())
    denominator = (merged.select(pl.max_horizontal("count", "count_2").sum()).item())
    if denominator == 0:
        return 0.0

    return numerator / denominator


def calc_dice_sorensen(merged):
    shared = ((pl.col("count") > 0) & (pl.col("count_2") > 0))
    intersection = (merged.filter(shared).height)
    sample1 = (merged.filter(pl.col("count") > 0).height)
    sample2 = (merged.filter(pl.col("count_2") > 0).height)
    denominator = sample1 + sample2
    if denominator == 0:
        return 0.0

    return (2 * intersection) / denominator



def compare_pair(pair):
    i, j, sample1_id, sample2_id = pair
    df1 = SAMPLES[sample1_id]
    df2 = SAMPLES[sample2_id]
    merged = compare_kmers(df1, df2)

    return (i, j, calc_jaccard(merged), calc_weighted_jaccard(merged), calc_dice_sorensen(merged))



def write_matrix(matrix, sample_ids, outfile):
    """
    Write only the upper triangle.
    The lower triangle is represented by NaN
    and will be filled by run_pipeline.py.
    """
    n = len(sample_ids)
    output = np.full((n, n), np.nan, dtype=np.float64)
    for i in range(n):
        for j in range(i, n):
            output[i, j] = matrix[i, j]

    df = pl.DataFrame(output, schema=sample_ids)
    df.insert_column(0, pl.Series("sample", sample_ids))
    tmpfile = outfile.with_name(outfile.name + ".tmp")
    df.write_csv(tmpfile, separator="\t")
    tmpfile.replace(outfile)



def process_batch(batch, jaccard_matrix, weighted_jaccard_matrix, dice_matrix, threads):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(compare_pair, pair) for pair in batch]
        for future in as_completed(futures):
            (i, j, jaccard, weighted_jaccard, dice) = future.result()
            jaccard_matrix[i, j] = jaccard
            weighted_jaccard_matrix[i, j] = weighted_jaccard
            dice_matrix[i, j] = dice

            yield 1


## Main
def main():
    parser = argparse.ArgumentParser(description=("Calculate pairwise k-mer similarities"))
    parser.add_argument("--indir", type=Path, required=True, help="Directory containing Parquet files")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Output prefix, e.g. 500k_27")
    parser.add_argument("--threads", type=int, default=16, help="Number of simultaneous comparisons")
    parser.add_argument("--min-count", type=int, default=5, help="Minimum k-mer count")
    parser.add_argument("--batch-size", type=int, default=10000, help=("Number of comparisons submitted at one time "))
    args = parser.parse_args()

    if args.threads < 1:
        parser.error("--threads must be >= 1")

    if args.batch_size < 1:
        parser.error("--batch-size must be >= 1")

    args.outdir.mkdir(parents=True, exist_ok=True)
    outputs = matrix_paths(args.outdir, args.prefix) if all(outfile.exists() for outfile in outputs):
        print(f"Pairwise comparison skipped, all matrices exists")

        return

    existing = [outfile.name for outfile in outputs if outfile.exists()]
    if existing:
        print(f"Pairwise comparison incomplete")
        for filename in existing:
            print(" Recalculating all matrices")
    else:
        print(f"{args.prefix}: matrices do not exist")

    parquet_files = sorted(args.indir.glob("*.parquet"))
    if not parquet_files:
        print(f"no Parquet files found in {args.indir}")
        raise SystemExit(1)

    print(f"Found {len(parquet_files)} Parquet files")
    SAMPLES.clear()

    for parquet_file in tqdm(parquet_files, desc="Loading k-mer database"):
        sample_id = parquet_file.stem
        SAMPLES[sample_id] = (pl.read_parquet(parquet_file).filter(pl.col("count") >= args.min_count))

    sample_ids = sorted(SAMPLES.keys())
    n_samples = len(sample_ids)
    print(f"Loaded {n_samples} samples")

    jaccard_matrix = np.zeros((n_samples, n_samples), dtype=np.float64)
    weighted_jaccard_matrix = np.zeros(n_samples, n_samples), dtype=np.float64)
    dice_matrix = np.zeros((n_samples, n_samples), dtype=np.float64)

    np.fill_diagonal(jaccard_matrix, 1.0)
    np.fill_diagonal(weighted_jaccard_matrix, 1.0)
    np.fill_diagonal(dice_matrix, 1.0)


    def pair_generator():
        for i in range(n_samples):
            for j in range(i + 1, n_samples):
                yield (i, j, sample_ids[i], sample_ids[j])

    total_comparisons = (n_samples * (n_samples - 1) // 2)
    print(f"{total_comparisons} total pairwise comparisons")
    print(f"Threads {args.threads}")
    print(f"Batch size {args.batch_size:,}")

    generator = pair_generator()
    processed = 0
    with tqdm(total=total_comparisons, desc="Calculating similarities") as pbar:
        while True:
            batch = []
            try:
                for _ in range(args.batch_size):
                    batch.append(next(generator))

            except StopIteration:
                pass

            if not batch:
                break

            for _ in process_batch(batch, jaccard_matrix, weighted_jaccard_matrix, dice_matrix, args.threads):
                processed += 1
                pbar.update(1)

    write_matrix(jaccard_matrix, sample_ids, outputs[0])
    write_matrix(weighted_jaccard_matrix, sample_ids, outputs[1])
    write_matrix(dice_matrix, sample_ids, outputs[2])
    for outfile in outputs:
        print(f"  {outfile}")

##
if __name__ == "__main__":
    main()
