#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script to computing pairwise snp distance betweeen isolates

Usage: python3 pairwise_snp_distance.py
'''

## import dependencies
import gzip
import glob
import os
import re
import csv
from collections import Counter

## set paths
BASE = "/lustre/BIF/nobackup/moret013/thesis_GEN"
VCF_DIR = os.path.join(BASE, "snelders_pipeline/missing_group_vcfs")
SUMMARY_IN = os.path.join(BASE, "clonal_group_summary.tsv")
SUMMARY_OUT = os.path.join(BASE, "clonal_group_pairwise_distances.tsv")


## define functions
def read_gt(sample_field, format_field):
    """Return a called genotype, or None if missing."""
    fmt = format_field.split(":")
    values = sample_field.split(":")
    if "GT" not in fmt:
        return None

    idx = fmt.index("GT")
    if idx >= len(values):
        return None

    gt = values[idx].replace("|", "/")
    if gt in ("", ".", "./.", ".|."):
        return None

    alleles = gt.split("/")
    if any(a == "." or a == "" for a in alleles):
        return None

    return tuple(sorted(alleles, key=int))



def analyze_vcf(path):
    total_snp_sites = 0
    total_discordant_pairs = 0
    total_comparable_pair_sites = 0
    sample_names = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                sample_names = line.rstrip("\n").split("\t")[9:]
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            ref = fields[3] ## count snps only
            alts = fields[4].split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue

            total_snp_sites += 1
            fmt = fields[8]
            genotypes = [read_gt(sample, fmt) for sample in fields[9:]]
            called = [gt for gt in genotypes if gt is not None]
            counts = Counter(called)
            n_called = len(called) ## find the number of pairs with snp calls at that site
            comparable_pairs = n_called * (n_called - 1) // 2
            identical_pairs = sum(n * (n - 1) // 2 for n in counts.values()) ## sum the pairs with identical genotypes
            discordant_pairs = comparable_pairs - identical_pairs
            total_comparable_pair_sites += comparable_pairs
            total_discordant_pairs += discordant_pairs

    n = len(sample_names)
    n_pairs = n * (n - 1) // 2

    if n_pairs == 0:
        raise ValueError(f"Fewer than two samples in {path}")

    observed_differences_per_pair = (total_discordant_pairs / n_pairs) ## observed diff. averaged over all possible pairs, missing pair/sites contribute no observed diff.
    difference_rate = (total_discordant_pairs / total_comparable_pair_sites if total_comparable_pair_sites else float("nan")) ##Fraction of the comparable pairs/sites calls that are not concordant
    standardized_distance = (difference_rate * total_snp_sites if total_comparable_pair_sites else float("nan")) ## standardize estimates and account for missingness using total observed diff. rate

    return {
        "group_size": n,
        "pair_count": n_pairs,
        "snp_sites": total_snp_sites,
        "discordant_pair_site_calls": total_discordant_pairs,
        "comparable_pair_site_calls": total_comparable_pair_sites,
        "observed_differences_per_pair": observed_differences_per_pair,
        "pairwise_snp_difference_rate": difference_rate,
        "standardized_mean_pairwise_distance": standardized_distance}

## Def main
def main():
    vcf_files = sorted(glob.glob(os.path.join(VCF_DIR, "ClonalGroup_*.vcf.gz")), key=lambda p: int(re.search(r"ClonalGroup_(\d+)", p).group(1)))
    if not vcf_files:
        raise FileNotFoundError(f"No group VCFs found in {VCF_DIR}")

    metadata = {}
    with open(SUMMARY_IN, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            metadata[row["clonal_group"]] = row

    output_rows = []

    for path in vcf_files:
        match = re.search(r"ClonalGroup_(\d+)", path)
        group = f"ClonalGroup_{match.group(1)}"
        result = analyze_vcf(path)
        meta = metadata.get(group, {})
        row = {"clonal_group": group, "MSH6_G233A_status": meta.get("MSH6_G233A_status", "Unknown"), "variant_count": meta.get("variant_count", "NA"), **result}
        output_rows.append(row)

    columns = [
        "clonal_group",
        "MSH6_G233A_status",
        "variant_count",
        "group_size",
        "pair_count",
        "snp_sites",
        "discordant_pair_site_calls",
        "comparable_pair_site_calls",
        "observed_differences_per_pair",
        "pairwise_snp_difference_rate",
        "standardized_mean_pairwise_distance"]

    with open(SUMMARY_OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Saved results at {SUMMARY_OUT}")

##
if __name__ == "__main__":
    main()

