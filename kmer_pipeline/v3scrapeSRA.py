#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 13/07/2026

Description:
    Python script for scraping sequence data from SRA/SRR numbers.
    Allows retrieval of an X amount of reads.

Usage:
    python3 scrapeSRA.py --sralist file.txt --num_reads int --outdir path
'''

## Import packages and API

import subprocess
from pathlib import Path
import requests
import random
import argparse
import tempfile
import gzip

ENA_api = "https://www.ebi.ac.uk/ena/portal/api/filereport"


## Define functions:

def get_runs(accession):
    """
    Return a list of run accessions.
    """
    prefixes = ("SRR", "ERR", "DRR")
    if accession.startswith(prefixes):
        return [accession]
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": "run_accession",
        "format": "tsv"}
    request = requests.get(ENA_api, params=params)
    request.raise_for_status()
    lines = request.text.strip().splitlines()
    if len(lines) < 2:
        return []

    return [line.split("\t")[0] for line in lines[1:]]



def fetch_fastq(accession):
    """
    Return FASTQ URLs from ENA.
    """
    parameters = {
        "accession": accession,
        "result": "read_run",
        "fields": "fastq_ftp",
        "format": "tsv"}
    request = requests.get(ENA_api, params=parameters)
    request.raise_for_status()
    lines = request.text.strip().splitlines()
    if len(lines) < 2:
        return []

    urls = lines[1].split("\t")[1]
    if not urls:
        return []

    return ["https://" + url if not url.startswith("http") else url for url in urls.split(";")]



def count_reads(fastq):
    """
    Count the number of reads in a gzipped FASTQ file.
    """
    count = 0
    with gzip.open(fastq, "rt") as f:
        for line_number, line in enumerate(f):
            if line_number % 4 == 0:
                count += 1
    return count



def subsample_paired_fastq(r1_file, r2_file, out_r1, out_r2, reads):
    """
    Randomly select paired reads while maintaining read pairing.
    """
    total_reads = count_reads(r1_file)
    if reads >= total_reads:
        print(f"Requested {reads} reads, but the file only contains {total_reads}, continuing with max reads")
        subprocess.run(["cp", str(r1_file), str(out_r1)], check=True)
        subprocess.run(["cp", str(r2_file), str(out_r2)], check=True)

        return

    selected_reads = set(random.sample(range(total_reads), reads))

    with gzip.open(r1_file, "rt") as r1, \
         gzip.open(r2_file, "rt") as r2, \
         gzip.open(out_r1, "wt") as out1, \
         gzip.open(out_r2, "wt") as out2:


        read_number = 0
        while True:
            rec1 = [r1.readline() for _ in range(4)]
            rec2 = [r2.readline() for _ in range(4)]

            if not rec1[0]:
                break

            if read_number in selected_reads:
                out1.writelines(rec1)
                out2.writelines(rec2)
            read_number += 1



def sample_fastq(urls, outfiles, reads=None):
    """
    Download FASTQ files and optionally subsample reads.

    Handles:
    - single-end data (1 FASTQ)
    - paired-end data (2 FASTQs)
    """

    seed = random.randint(1, 100000)
    random.seed(seed)
    if len(urls) == 1:
        print("Detected single-end FASTQ")
        if reads is None:
            cmd = f"""curl -L {urls[0]} | gzip -c > {outfiles[0]}"""
        else:
            cmd = f"""curl -L {urls[0]} | gunzip | seqtk sample -s{seed} - {reads} | gzip > {outfiles[0]}"""

        subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)

    elif len(urls) == 2:
        print("Detected paired-end FASTQ")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_r1 = Path(tmp) / "R1.fastq.gz"
            tmp_r2 = Path(tmp) / "R2.fastq.gz"
            subprocess.run(f"curl -L {urls[0]} -o {tmp_r1}", shell=True, check=True)
            subprocess.run(f"curl -L {urls[1]} -o {tmp_r2}", shell=True, check=True)
            if reads is None:
                subprocess.run(["cp", str(tmp_r1), str(outfiles[0])], check=True)
                subprocess.run(["cp", str(tmp_r2), str(outfiles[1])], check=True)
            else:
                subsample_paired_fastq(tmp_r1, tmp_r2, outfiles[0], outfiles[1], reads)

    else:
        raise ValueError("Unexpected number of FASTQ files returned by ENA")


## Main

def main():
    parser = argparse.ArgumentParser(description="Fetch and subsample ENA/SRA/SRR FASTQ files")
    parser.add_argument("--sralist", required=True, help="Text file containing SRA accessions")
    parser.add_argument("--num_reads", type=int, default=None, help="Number of reads to retain (default= all reads)")
    parser.add_argument("--outdir", default="Fetched_reads", help="Output directory")
    args = parser.parse_args()

    if args.num_reads is not None and args.num_reads <= 0:
        raise ValueError("Number of reads must be greater than zero")

    outdir = Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    with open(args.sralist) as f:
        accessions = [x.strip() for x in f if x.strip()]

    for accession in accessions:
        runs = get_runs(accession)
        for run in runs:
            try:
                urls = fetch_fastq(run)
                if not urls:
                    print(f"No FASTQ found for {run}")
                    continue

                if len(urls) == 2:
                    outfiles = [
                        outdir / f"{run}_R1.fastq.gz",
                        outdir / f"{run}_R2.fastq.gz"]
                else:
                    outfiles = [outdir / f"{run}.fastq.gz"]

                sample_fastq(urls, outfiles, args.num_reads)
                print(f"Finished {run}")

            except subprocess.CalledProcessError:
                print(f"FAILED processing {run}")
                continue

            except requests.exceptions.RequestException:
                print(f"FAILED querying ENA: {run}")
                continue


if __name__ == "__main__":
    main()
