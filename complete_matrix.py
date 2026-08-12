#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 30/07/2026

Description:
    Mirrors the produced triangular matrices produced by pairwise_comparison.py
    since only the upper triangle is computed, via mirroring, we save time

Usage: python3 jellyfish_to_parquet.py --indir /path/to/files.jf --outdir /path/to/outdir
'''
## Import statements
from pathlib import Path
import polars as pl
import numpy as np
from tqdm import tqdm

## Hardcoded variables
indir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/kmerdatabase/outfiles")
outdir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/kmerdatabase/matrices_full")
outdir.mkdir(parents=True,exist_ok=True)

## Define functions

def mirror_matrix(infile, outfile):
    df = pl.read_csv(infile, separator ="\t")
    sample_ids = df["sample"].to_list()
    matrix = df.drop("sample").to_numpy()
    matrix = np.where(np.isnan(matrix), matrix.T, matrix)

    out_df = pl.DataFrame(matrix, schema=sample_ids)
    out_df.insert_column(0, pl.Series("sample", sample_ids))
    out_df.write_csv(outfile, separator = "\t")

## Main
def main():
    tsv_files = list(indir.glob("*.tsv"))
    for infile in tqdm(tsv_files, desc= "Mirroring matrices"):
        outfile = outdir/infile.name
        mirror_matrix(infile, outfile)

##
if __name__ == "__main__":
    main()
