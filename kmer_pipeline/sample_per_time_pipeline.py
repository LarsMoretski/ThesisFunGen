#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026
Description:
    Python script for automating the kmer pipeline

Usage: python3 sample_per_time_pipeline.py
'''

## import packages
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

## Set paths
BASE_DIR = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/final_analysis")
RAW_READS = BASE_DIR / "01_Raw_Reads"
SUBSAMPLING = BASE_DIR / "02_Subsampling"
JELLYFISH = BASE_DIR / "03_Jellyfish"
PARQUET = BASE_DIR / "04_Parquet_files"
RESULTS_DIR = BASE_DIR / "06_Results"
SAMPLE_LIST = BASE_DIR / "00_raw_input" / "final_isolate_list.txt"
SUBSAMPLER = Path("/lustre/BIF/nobackup/moret013/thesis_GENfinalscripts/finalrun/subsampler.py")
JELLYFISH_SCRIPT = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/finalscripts/finalrun/kmercounting.sh")
CONVERTER = Path("/lustre/BIF/nobackup/moret013/thesis_GEN/finalscripts/finalrun/jellyfish_to_parquet.py")

## Set read sizes and value for K
READ_COUNTS = [2_000_000]
KMER_SIZES = [27]

## Set seed, number of threads and the minimum count for the kmers
SEED = 42
SUBSAMPLER_THREADS = 1
JELLYFISH_THREADS = 60
MIN_COUNT = 5
FAILED_LOG = RESULTS_DIR / "failed_samples.txt"
MISSING_RAW_LOG = RESULTS_DIR / "missing_raw_pipeline.log"


## Define functions
for directory in [SUBSAMPLING, JELLYFISH, PARQUET, RESULTS_DIR,]:
    directory.mkdir(parents=True, exist_ok=True)



def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")



def log_failed(sample, reads, kmer, step, message):
    with open(FAILED_LOG, "a") as f:
        f.write(
            f"{timestamp()}\t"
            f"{sample}\t"
            f"{reads}\t"
            f"{kmer}\t"
            f"{step}\t"
            f"{message}\n")


def log_missing_raw(sample, reads, kmer, missing_r1, missing_r2):
    with open(MISSING_RAW_LOG, "a") as f:
        f.write(
            f"{timestamp()}\t"
            f"{sample}\t"
            f"{reads}\t"
            f"{kmer}\t"
            f"R1={'MISSING' if missing_r1 else 'OK'}\t"
            f"R2={'MISSING' if missing_r2 else 'OK'}\n")



def read_sample_list():
    if not SAMPLE_LIST.exists():
        raise FileNotFoundError(f"Sample list not found")

    samples = []
    with open(SAMPLE_LIST) as f:
        for line in f:
            sample = line.strip()

            if not sample:
                continue

            samples.append(sample)

    return samples



def raw_reads_exist(sample):
    r1 = RAW_READS / f"{sample}_R1.fastq.gz"
    r2 = RAW_READS / f"{sample}_R2.fastq.gz"
    r1_exists = r1.exists() and r1.stat().st_size > 0
    r2_exists = r2.exists() and r2.stat().st_size > 0

    return r1_exists, r2_exists



def get_subsample_paths(sample, reads):
    r1 = SUBSAMPLING / f"{sample}_R1.fastq.gz"
    r2 = SUBSAMPLING / f"{sample}_R2.fastq.gz"

    return r1, r2


def get_jellyfish_path(sample, reads, kmer):
    directory = JELLYFISH / f"{reads}_reads" / f"k{kmer}"
    return directory / f"{sample}.jf"



def get_parquet_path(sample, reads, kmer):
    directory = PARQUET / f"{reads}_reads" / f"k{kmer}"
    return directory / f"{sample}.parquet"


def run_command(command, log_file=None):
    if log_file is not None:
        with open(log_file, "a") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    else:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    return result


def process_configuration(sample, reads, kmer):
    subsample_r1, subsample_r2 = get_subsample_paths(sample, reads)
    jellyfish_file = get_jellyfish_path(sample, reads, kmer)
    parquet_file = get_parquet_path(sample, reads, kmer)

    if parquet_file.exists() and parquet_file.stat().st_size > 0:
        print(f"parquet already exists:")

        return "SUCCESS"


    if subsample_r1.exists() and subsample_r2.exists():
        print("subsampled reads already exist.")

    else:
        subsample_r1.parent.mkdir(parents=True, exist_ok=True,)
        command = [sys.executable, str(SUBSAMPLER), "--indir", str(RAW_READS), "--outdir", str(SUBSAMPLING), "--reads", str(reads), "--threads", str(SUBSAMPLER_THREADS),
            "--seed", str(SEED), "--sample", sample, "--sralist", str(SAMPLE_LIST), "--parquetdir", str(PARQUET)]

        result = run_command(command)
        if result.returncode != 0:
            if result.stdout:
                print(result.stdout)

            log_failed(sample, reads, kmer, "STEP1_SUBSAMPLING","subsampling command failed")
            if subsample_r1.exists():
                subsample_r1.unlink()

            if subsample_r2.exists():
                subsample_r2.unlink()

            return "FAILED"

        if not subsample_r1.exists() or not subsample_r2.exists():
            print("[FAILED] Subsampling completed but output is missing.")

            log_failed(sample, reads, kmer, "STEP1_SUBSAMPLING", "subsampled FASTQ files were not created")
            return "FAILED"

        print("subsampling complete.")

## jellyfish
    if jellyfish_file.exists() and jellyfish_file.stat().st_size > 0:
        print("Jellyfish file already exists.")
    else:
        jellyfish_file.parent.mkdir(parents=True, exist_ok=True,)
        command = ["bash", str(JELLYFISH_SCRIPT), "--indir", str(subsample_r1.parent), "--outdir", str(jellyfish_file.parent), "--sample", sample, 
        "--k", str(kmer), "--threads", str(JELLYFISH_THREADS)]

        result = run_command(command)
        if result.returncode != 0:
            print("[FAILED] Jellyfish counting.")
            if result.stdout:
                print(result.stdout)

            log_failed(sample, reads, kmer, "STEP2_JELLYFISH", "Jellyfish command failed")
            return "FAILED"

        if not jellyfish_file.exists() or jellyfish_file.stat().st_size == 0:
            print("Jellyfish completed but output is missing.")
            log_failed(sample, reads, kmer, "STEP2_JELLYFISH", "Jellyfish file was not created")

            return "FAILED"

        print("Jellyfish counting complete.")

## Jellyfish to parquet
    parquet_file.parent.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(CONVERTER), "--indir", str(jellyfish_file.parent), "--outdir", str(parquet_file.parent), "--sample", sample,]
    result = run_command(command)
    if result.returncode != 0:
        print("Jellyfish to Parquet conversion.")

        if result.stdout:
            print(result.stdout)

        log_failed(sample, reads, kmer, "STEP3_PARQUET", "parquet conversion failed")
        return "FAILED"

    if not parquet_file.exists() or parquet_file.stat().st_size == 0:
        log_failed(sample, reads, kmer, "STEP3_PARQUET", "Expected Parquet file was not created")
        return "FAILED"

    print("Parquet conversion complete.")
    print("Removing temporary files")
    if jellyfish_file.exists():
        jellyfish_file.unlink()

    if subsample_r1.exists():
        subsample_r1.unlink()

    if subsample_r2.exists():
        subsample_r2.unlink()

    print("Sample complete")

    return "SUCCESS"


## Def main
def main():
    samples = read_sample_list()
    successful = 0
    missing_raw = 0
    failed = 0
    for sample_number, sample in enumerate(samples, start=1):
        sample_success = True
        for reads in READ_COUNTS:
            for kmer in KMER_SIZES:
                status = process_configuration(sample, reads, kmer)
                if status == "MISSING_RAW":
                    missing_raw += 1
                    sample_success = False
                    break

                elif status == "FAILED":
                    failed += 1
                    sample_success = False
                    break

                elif status == "SUCCESS":
                    successful += 1

            if not sample_success:
                break

    print(f"Successful configurations: {successful}")
    print(f"Missing raw reads:         {missing_raw}")
    print(f"Failed configurations:     {failed}")
    print(f"Failed log at {FAILED_LOG}")
    print(f"Missing raw-read log at {MISSING_RAW_LOG}")


##
if __name__ == "__main__":
    main()
