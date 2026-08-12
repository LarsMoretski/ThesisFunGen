#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/07/2026

Description:
    Python script for parsing jellyfish output files (.jf) to .parquet format.
    Generates the .parquet files used by DuckDB

Usage: python3 jellyfish_to_parquet.py --indir /path/to/files.jf --outdir /path/to/outdir
'''

## Import packages
import polars as pl
import argparse
import subprocess
from pathlib import Path

## Define functions
def convert_to_parquet(jf_file, outdir):
    sample=jf_file.stem
    outfile= outdir/f"{sample}.parquet"
    if outfile.exists():
        print(f"{sample}.parquet already exists, skipping {sample}...")
        return

    result = subprocess.run(["jellyfish", "dump", "-c", str(jf_file)], capture_output=True, text=True, check=True)
    rows = [line.split() for line in result.stdout.splitlines()]
    df = pl.DataFrame(rows, schema=["kmer", "count"])
    df = df.with_columns(pl.col("count").cast(pl.UInt32))
    df.write_parquet(outfile)
    print(f" Created {outfile}")

## Main
def main():
    indir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/subdataset/subset_unk_tryout")
    outdir = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/kmerdatabase/subset_fullreads")
    outdir.mkdir(parents=True, exist_ok=True)

    for jf_file in indir.glob("*.jf"):
        convert_to_parquet(jf_file, outdir)


if __name__ == "__main__":
    main()
