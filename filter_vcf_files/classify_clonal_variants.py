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
import csv
import gzip
from pathlib import Path
from collections import Counter

# Input files
GENOTYPE_FILE = Path("msh6_G233A_genotypes.tsv")
VCF_DIR = Path("snelders_pipeline/missing_group_vcfs")
OUTPUT_FILE = Path("clonal_vcf_MSH6_G233A_summary.tsv")

# Read isolate-level MSH6 G233A genotypes
genotypes = {}
with GENOTYPE_FILE.open("r") as handle:
    reader = csv.reader(handle, delimiter="\t")
    for row in reader:
        if len(row) < 2:
            continue

        sample = row[0].strip()
        gt = row[1].strip()

        if sample.lower() in {"sample", "sample_id", "isolate", "isolate_id"}:
            continue

        if gt == "1":
            genotypes[sample] = "Present"
        elif gt == "0":
            genotypes[sample] = "Absent"
        else:
            genotypes[sample] = "Unknown"

## Define functions
def get_vcf_samples(vcf_path):
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                return columns[9:] if len(columns) > 9 else []

    raise ValueError(f"No #CHROM header found in {vcf_path}")


vcf_files = sorted(VCF_DIR.glob("ClonalGroup_*.vcf.gz"))
if not vcf_files:
    raise FileNotFoundError(f"No ClonalGroup_*.vcf.gz files found in {VCF_DIR}")

results = []
for i, vcf in enumerate(vcf_files, start=1):
    samples = get_vcf_samples(vcf)
    statuses = [genotypes.get(sample, "Unknown") for sample in samples]
    counts = Counter(statuses)
    present = counts["Present"]
    absent = counts["Absent"]
    unknown = counts["Unknown"]

    if present > 0 and absent > 0:
        group_status = "Mixed"
    elif present > 0:
        group_status = "Present"
    elif absent > 0:
        group_status = "Absent"
    else:
        group_status = "Unknown"

    results.append({
        "clonal_group": vcf.stem.removesuffix(".vcf"),
        "sample_count": len(samples),
        "G233A_present": present,
        "G233A_absent": absent,
        "unknown_or_missing": unknown,
        "group_status": group_status,
        "unmatched_samples": sum(
            sample not in genotypes for sample in samples),})

    if i % 25 == 0:
        print(f"Processed {i}/{len(vcf_files)} clonal groups...")

if results:
    with OUTPUT_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=results[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


