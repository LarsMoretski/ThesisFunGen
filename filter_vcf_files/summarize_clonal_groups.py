#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script for summarizing the clonal groups


Usage: python3 summarize_clonal_groups.py
'''

## import packages
import gzip
import glob
import os
import re
import csv

## set paths
BASE_DIR = "/lustre/BIF/nobackup/moret013/thesis_GEN"
GROUP_VCF_DIR = os.path.join(BASE_DIR, "snelders_pipeline")
VARIABLE_VCF_DIR = os.path.join(BASE_DIR, "clonal_group_variable_vcfs")
GENOTYPE_FILE = os.path.join(BASE_DIR, "msh6_G233A_genotypes.tsv")
OUTPUT_FILE = os.path.join(BASE_DIR, "clonal_group_summary.tsv")


## define functions
def open_vcf(path):
    return gzip.open(path, "rt")



def get_samples_and_size(path):
    """Get isolate names from a clonal group's VCF header."""
    with open_vcf(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                columns = line.rstrip("\n").split("\t")
                samples = columns[9:]
                return samples, len(samples)

    raise ValueError(f"No #CHROM header found in {path}")



def count_variants(path):
    """Count variant records, excluding VCF headers."""
    count = 0
    with open_vcf(path) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                count += 1

    return count



def load_msh6_genotypes(path):
    """Read sample ID and genotype (0, 1, or missing)."""
    genotypes = {}
    with open(path) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 2:
                continue

            sample, gt = fields[0], fields[1]

            if gt not in {"0", "1", "."}:
                continue

            genotypes[sample] = gt

    return genotypes



def classify_group(samples, genotypes):
    calls = [genotypes.get(sample, ".") for sample in samples]
    called = [gt for gt in calls if gt in {"0", "1"}]
    if not called:
        return "Unknown"

    if "0" in called and "1" in called:
        return "Mixed"

    if all(gt == "1" for gt in called):
        return "A233"

    return "G233"



def main():
    genotypes = load_msh6_genotypes(GENOTYPE_FILE)
    group_files = sorted(glob.glob(os.path.join(GROUP_VCF_DIR, "ClonalGroup_*.vcf.gz")), key=lambda p: int(re.search(r"ClonalGroup_(\d+)", p).group(1)))
    if not group_files:
        raise FileNotFoundError(f"No clonal group VCFs found in {GROUP_VCF_DIR}")

    rows = []
    for group_path in group_files:
        match = re.search(r"ClonalGroup_(\d+)", group_path)
        group_id = f"ClonalGroup_{match.group(1)}"
        samples, group_size = get_samples_and_size(group_path)
        variable_path = os.path.join(VARIABLE_VCF_DIR, f"{group_id}.variable.vcf.gz")

        if os.path.exists(variable_path):
            variant_count = count_variants(variable_path)
        else:
            # Missing file is not the same as zero variants, so they get marked differently
            variant_count = "NA"

        status = classify_group(samples, genotypes)

        rows.append({
            "clonal_group": group_id,
            "variant_count": variant_count,
            "MSH6_G233A_status": status,
            "group_size": group_size})

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
                "clonal_group",
                "variant_count",
                "MSH6_G233A_status",
                "group_size",], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Processed {len(rows)} clonal groups.")
    print(f"Output at {OUTPUT_FILE}")

##
if __name__ == "__main__":
    main()

