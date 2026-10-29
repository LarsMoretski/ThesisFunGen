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
from pathlib import Path
import csv
from collections import Counter, defaultdict

## Set directories
base = Path("mutational_patterns")
matrix = base / "combined/All_isolates.SBS96.tsv"  ## this is where the SBS96 matrices of all isolates go
genotype_file = Path("msh6_G233A_genotypes.tsv")   ## this is the genotyping file containing the status of pos 233
vcf_root = base / "all_group_vcfs"                 ## this is the dir where all separated VCF files are present
outdir = base / "combined"
outdir.mkdir(parents=True, exist_ok=True)

## load the genotypes
genotypes = {}
with genotype_file.open() as f:
    for line in f:
        fields = line.split()
        if len(fields) >= 2:
            genotypes[fields[0]] = fields[1]

# Read isolate names from the SBS96 matrix
with matrix.open() as f:
    samples = f.readline().rstrip("\n").split("\t")[1:]

if len(samples) != len(set(samples)):
    raise SystemExit("ERROR: duplicate sample names in SBS96 matrix")

sample_groups = defaultdict(list)
for group_dir in sorted(vcf_root.glob("ClonalGroup_*")):
    if not group_dir.is_dir():
        continue

    for vcf in group_dir.glob("*.vcf.gz"):
        sample = vcf.name.removesuffix(".vcf.gz")
        sample_groups[sample].append(group_dir.name)


missing_group = [s for s in samples if s not in sample_groups]
multiple_groups = {s: sample_groups[s] for s in samples if len(sample_groups[s]) > 1}
missing_genotype = [s for s in samples if s not in genotypes]

if missing_group or multiple_groups or missing_genotype:
    print("Samples without group:", len(missing_group), missing_group[:10])
    print("Samples in multiple groups:", len(multiple_groups), list(multiple_groups.items())[:10])
    print("Samples without genotype entry:", len(missing_genotype), missing_genotype[:10])
    raise SystemExit("ERROR: resolve mapping issues before analysis")

def label(g):
    return {"0": "G233", "1": "A233", ".": "Unknown"}.get(g, "Unknown")

outfile = outdir / "All_isolates_metadata.tsv"
with outfile.open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["sample_id", "clonal_group", "genotype_code", "MSH6_genotype"])
    for sample in samples:
        group = sample_groups[sample][0]
        g = genotypes[sample]
        writer.writerow([sample, group, g, label(g)])

print("Wrote:", outfile)
print("Total isolates:", len(samples))
print("Genotype counts:", dict(Counter(label(genotypes[s]) for s in samples)))
print("Clonal groups represented:", len({sample_groups[s][0] for s in samples}))
print("Group-size distribution (isolates per group):")
group_counts = Counter(sample_groups[s][0] for s in samples)
print(dict(sorted(Counter(group_counts.values()).items())))
