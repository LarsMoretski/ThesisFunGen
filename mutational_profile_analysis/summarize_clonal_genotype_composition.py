#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    python script to

Usage: python3 jellyfish_to_parquet.py --indir /path/to/files.jf --outdir /path/to/outdir
'''

## import packages
import pandas as pd

## load data and set path
metadata_file = "mutational_patterns/combined/All_isolates_metadata.tsv"
outdir = "mutational_patterns/combined"


meta = pd.read_csv(metadata_file, sep="\t")
rows = []
for group, df in meta.groupby("clonal_group"):
    counts = df["MSH6_genotype"].value_counts()
    n_a = int(counts.get("A233", 0))
    n_g = int(counts.get("G233", 0))
    n_unknown = int(counts.get("Unknown", 0))

    if n_a > 0 and n_g > 0:
        status = "Both A233 and G233"
    elif n_a > 0:
        status = "A233 only"
    elif n_g > 0:
        status = "G233 only"
    else:
        status = "Unknown only"

    rows.append({
        "clonal_group": group,
        "group_size": len(df),
        "A233_count": n_a,
        "G233_count": n_g,
        "unknown_count": n_unknown,
        "genotype_composition": status})

result = pd.DataFrame(rows).sort_values("clonal_group")
result.to_csv(f"{outdir}/ClonalGroup_genotype_composition.tsv", sep="\t", index=False)

print("Number of groups by genotype composition")
print(result["genotype_composition"].value_counts().to_string())
print("Isolates represented in each group category")
print(result.groupby("genotype_composition")["group_size"].agg(["count", "sum", "median", "min", "max"]).to_string())
print("Groups containing both genotypes")
mixed = result[result["genotype_composition"] == "Both A233 and G233"]
print(mixed.to_string(index=False) if len(mixed) else "None")
print(f"Saved at {outdir}/ClonalGroup_genotype_composition.tsv")
