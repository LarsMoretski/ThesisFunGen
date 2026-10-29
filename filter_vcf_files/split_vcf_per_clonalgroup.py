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

#import packages
import argparse
import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

## Define functions
def run(command, capture=False):
    result = subprocess.run(command, check=True, text=True, stdout=subprocess.PIPE if capture else None,)
    return result.stdout if capture else None



def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


## Def main
def main():
    parser = argparse.ArgumentParser(description="Split multi-sample VCF into clonal-group VCFs")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ")
    parser.add_argument("--isolates", required=True, help="Clonal isolate TSV")
    parser.add_argument("--outdir", default="clonal_group_vcfs")
    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--group-col", default="clonal_group")
    parser.add_argument("--strip-prefix", action="store_true", help="Also try matching the final underscore-delimited ID")
    args = parser.parse_args()

    vcf = Path(args.vcf).resolve()
    isolate_file = Path(args.isolates).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    with open(isolate_file, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            sys.exit("ERROR: The isolate TSV has no header")

        for col in (args.sample_col, args.group_col):
            if col not in reader.fieldnames:
                sys.exit(f"ERROR: Column '{col}' not found ")

        records = []
        for row in reader:
            sample = (row[args.sample_col] or "").strip()
            group = (row[args.group_col] or "").strip()
            if sample and group:
                records.append((sample, group))

    if not records:
        sys.exit("ERROR: No sample/group records found.")

    vcf_samples = run(["bcftools", "query", "-l", str(vcf)], capture=True).splitlines()
    vcf_sample_set = set(vcf_samples)

    def match_sample(sample):
        if sample in vcf_sample_set:
            return sample
        if args.strip_prefix:
            suffix = sample.rsplit("_", 1)[-1]
            if suffix in vcf_sample_set:
                return suffix
        return None

    matched = []
    unmatched = []
    seen = set()
    for sample, group in records:
        vcf_id = match_sample(sample)
        if vcf_id is None:
            unmatched.append((sample, group))
            continue

        key = (vcf_id, group)
        if key not in seen:
            matched.append((vcf_id, sample, group))
            seen.add(key)

    unmatched_file = outdir / "unmatched_samples.tsv"
    with open(unmatched_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id_from_table", "clonal_group"])
        writer.writerows(sorted(unmatched, key=lambda x: (x[1], x[0])))

    if not matched:
        sys.exit(f"ERROR: No isolate IDs matched the VCF")

    sorted_table = outdir / "clonal_isolates_sorted.tsv"
    with open(sorted_table, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "vcf_sample_id", "clonal_group"])
        for vcf_id, original_id, group in sorted(matched, key=lambda x: (x[2], x[1])):
            writer.writerow([original_id, vcf_id, group])

    groups = {}
    for vcf_id, original_id, group in matched:
        groups.setdefault(group, set()).add(vcf_id)

    summary_file = outdir / "clonal_group_summary.tsv"
    with open(summary_file, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["clonal_group", "isolates_matched", "vcf_file"])

        for group in sorted(groups):
            samples = sorted(groups[group])
            prefix = safe_name(group)
            sample_file = outdir / f".{prefix}.samples.txt"
            output_vcf = outdir / f"{prefix}.vcf.gz"
            sample_file.write_text("\n".join(samples) + "\n", encoding="utf-8")
            if output_vcf.exists():
                sys.exit(f"ERROR: {output_vcf} already exists")

            run(["bcftools", "view", "-S", str(sample_file), "-Oz", "-o", str(output_vcf), str(vcf)])
            run(["bcftools", "index", "-t", str(output_vcf)])
            sample_file.unlink()

            writer.writerow([group, len(samples), output_vcf.name])
            print(f"{group}: {len(samples)} isolates -> {output_vcf.name}")

##
if __name__ == "__main__":
    main()

