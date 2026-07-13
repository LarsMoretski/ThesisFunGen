#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 13/07/2026
Description: Python script for scraping sequence data from SRA/SRR numbers
             This script allows for an XX percent of the data to be retrieved
Usage: python3 scrapeSRA.py --sralist file.txt --percentage int
'''

## Import packages and statements:

import subprocess
from pathlib import Path
import requests
import random
import argparse

ENA_api = "https://www.ebi.ac.uk/ena/portal/api/filereport"

## Define functions:


def get_runs(accession):
    """
    Description: return a list of run accessions from a project accession.
    """

    prefixes = ("SRR", "ERR", "DRR")

    if accession.startswith(prefixes):
        return [accession]

    params = {"accession": accession,
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
    Description: Fetches the URLs of the run accessions, queries the ENA api, and lists the fastq name.  
    """
    parameters = {"accession": accession,
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

    return [
        "https://" + url if not url.startswith("http") else url
        for url in urls.split(";")
    ]



def sample_fastq(url, outfile, fraction):
    """
    Description: downloads the URL, and takes an x percentage of the file using seqtk
                 with a predefined percentage/fraction. The intermediate files will be
                 removed, and only the fraction file is kept.
    """
    seed = random.randint(1, 100000)

    cmd = f"""
    curl -L {url} |
    gunzip |
    seqtk sample -s{seed} - {fraction} |
    gzip > {outfile}
    """

    subprocess.run(cmd, shell=True, check=True)




def concatenate(files, outfile):
    """
    Description: concatenates Illumina fwd (r1), and rv (r2) reads into one file.
                 deletes the unmerged one.
    """
    with open(outfile, "wb") as out:
         subprocess.run(["cat"] + [str(x) for x in files], stdout=out, check=True)



## Main
def main():
    parser = argparse.ArgumentParser(description = "Fetch and subsample ENA/SRA/SRR fastq files")
    parser.add_argument("--sralist", required=True, help="Text file containing the SRA accessions")
    parser.add_argument("--percentage", type=float, default = 100, help="Percentage of the reads take from the files (default = 100%)")
    parser.add_argument("--outdir", default = "Fetched_reads", help= "Output directory of the fetched files")
    args = parser.parse_args()

    fraction = args.percentage / 100

    outdir= Path(args.outdir)
    outdir.mkdir(exist_ok=True)

    if args.percentage <= 0 or args.percentage > 100:
        raise ValueError("Percentage must be between 0 and 100")

    with open(args.sralist) as f:
        accessions = [x.strip() for x in f if x.strip()]

    for accession in accessions:
        runs = get_runs(accession)
        for run in runs:
            try:
                urls = fetch_fastq(run)
                if not urls:
                    continue

                fetched_files = []

                for i, url in enumerate(urls):
                    outfile = outdir / f"{run}_{i+1}.subset.fastq.gz"
                    sample_fastq(url, outfile, fraction)
                    fetched_files.append(outfile)

                if len(fetched_files) == 2:
                    merged = outdir / f"{run}.fastq.gz"
                    concatenate(fetched_files, merged)
                    for file in fetched_files:
                        file.unlink()

                elif len(fetched_files) == 1:
                    fetched_files[0].rename(outdir / f"{run}.fastq.gz")

            except subprocess.CalledProcessError:
                print(f"FAILED during download/subsampling: {run}")
                continue

            except requests.exceptions.RequestException:
                print(f"FAILED querying ENA: {run}")
                continue



if __name__ == "__main__":
    main()
