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
import gzip
import glob
import os
import re
import csv
from collections import Counter

## Set paths
BASE = "/lustre/BIF/nobackup/moret013/thesis_GEN"
VCF_DIR = os.path.join(BASE, "snelders_pipeline/missing_group_vcfs")
SUMMARY_IN = os.path.join(BASE, "clonal_group_summary.tsv")
SUMMARY_OUT = os.path.join(BASE, "clonal_group_pairwise_distances_untrimmed.tsv")


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
    if any(a in ("", ".") for a in alleles):
        return None

    try:
        return tuple(sorted(alleles, key=int))
    except ValueError:
        return None


def analyze_vcf(path):
    total_snp_sites = 0
    total_discordant = 0
    total_comparable = 0
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
            ref = fields[3]
            alts = fields[4].split(",")

            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                continue

            total_snp_sites += 1
            fmt = fields[8]
            genotypes = [read_gt(sample, fmt) for sample in fields[9:]]
            called = [gt for gt in genotypes if gt is not None]
            counts = Counter(called)
            n_called = len(called)
            comparable_pairs = n_called * (n_called - 1) // 2
            identical_pairs = sum(n * (n - 1) // 2 for n in counts.values())
            total_comparable += comparable_pairs
            total_discordant += comparable_pairs - identical_pairs

    n = len(sample_names)
    n_pairs = n * (n - 1) // 2

    if n_pairs == 0:
        raise ValueError(f"Fewer than two samples in {path}")

    possible_pair_site_comparisons = n_pairs * total_snp_sites
    difference_rate = (total_discordant / total_comparable if total_comparable else float("nan")
    callability = (total_comparable / possible_pair_site_comparisons if possible_pair_site_comparisons else float("nan"))

    return {
        "group_size": n,
        "pair_count": n_pairs,
        "snp_sites": total_snp_sites,
        "discordant_pair_site_calls": total_discordant,
        "comparable_pair_site_calls": total_comparable,
        "observed_differences_per_pair": total_discordant / n_pairs,
        "pairwise_snp_difference_rate": difference_rate,
        "callability": callability,
        "standardized_mean_pairwise_distance": (
            difference_rate * total_snp_sites
            if total_comparable else float("nan"),}

## Def main
def main():
    vcf_files = sorted(
        glob.glob(os.path.join(VCF_DIR, "ClonalGroup_*.vcf.gz")),
        key=lambda p: int(re.search(r"ClonalGroup_(\d+)", p).group(1)))

    if not vcf_files:
        raise FileNotFoundError(f"No VCFs found in {VCF_DIR}")

    metadata = {}
    with open(SUMMARY_IN, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            metadata[row["clonal_group"]] = row

    rows = []
    seen = set()

    for i, path in enumerate(vcf_files, start=1):
        match = re.search(r"ClonalGroup_(\d+)", path)
        group = f"ClonalGroup_{int(match.group(1)):03d}"

        if group in seen:
            raise ValueError(f"Duplicate group ID: {group}")
        seen.add(group)

        result = analyze_vcf(path)
        meta = metadata.get(group, {})

        rows.append({
            "clonal_group": group,
            "MSH6_G233A_status": meta.get("MSH6_G233A_status", "Unknown"),
            "variant_count": meta.get("variant_count", "NA"), **result,})

        print(
            f"[{i}/{len(vcf_files)}] {group}: "
            f"n={result['group_size']}, "
            f"observed differences/pair="
            f"{result['observed_differences_per_pair']:.3f}, "
            f"rate={result['pairwise_snp_difference_rate']:.8f}, "
            f"callability={result['callability']:.4f}")

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
        "callability",
        "standardized_mean_pairwise_distance",]

    with open(SUMMARY_OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

