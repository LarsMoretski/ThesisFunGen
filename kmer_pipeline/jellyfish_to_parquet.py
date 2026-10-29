#!/usr/bin/env python3

"""
Author: Lars Moret, MSc student @ FunGen

Description:
    Convert ONE Jellyfish database (.jf) to a Parquet file.

Resume behavior:
    - Existing .parquet -> skip
    - Missing .parquet -> convert
    - Failed conversion -> retry on next run
"""

## import packages
import argparse
import subprocess
from pathlib import Path
import polars as pl


##Define functions
def convert_to_parquet(jf_file, outdir):
    sample = jf_file.stem
    outfile = (outdir / f"{sample}.parquet")
    if outfile.exists():
        print(f"{outfile.name} already exists")

        return

    tmpfile = outfile.with_name(outfile.name + ".tmp")

    try:
        print(f"Converting {jf_file.name}")
        result = subprocess.run(["jellyfish", "dump", "-c", str(jf_file)], capture_output=True, text=True, check=True)
        rows = [line.split() for line in result.stdout.splitlines() if line.strip()]
        if not rows:
            raise ValueError(f"No k-mers found in {jf_file}")

        df = pl.DataFrame(rows, schema=["kmer", "count"], orient="row")
        df = df.with_columns(pl.col("count").cast(pl.UInt32))
        df.write_parquet(tmpfile)

        if not tmpfile.exists():
            raise RuntimeError("Parquet file was not created")

        if tmpfile.stat().st_size == 0:
            raise RuntimeError("Parquet file is empty")

        tmpfile.replace(outfile

    except Exception:
        if tmpfile.exists():
            tmpfile.unlink()

        raise


## Def main
def main():
    parser = argparse.ArgumentParser(description=("Convert one Jellyfish file to Parquet"))
    parser.add_argument("--indir", type=Path, required=True, help="dir containing .jf file")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    parser.add_argument("--sample", required=True, help="Sample name")
    args = parser.parse_args()

    jf_file = (args.indir / f"{args.sample}.jf")
    if not jf_file.exists():
        raise FileNotFoundError(f"{jf_file} file does not exist")

    args.outdir.mkdir(parents=True, exist_ok=True)
    convert_to_parquet(jf_file, args.outdir)


##
if __name__ == "__main__":
    main()
