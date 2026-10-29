#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script for summarizing the SBS96 counts

Usage: python3 summarize_sbs96_counts.py
'''

## import packages
import pandas as pd


## load the metadata and the SBS96 mutational profile matrix
matrix_file = "mutational_patterns/combined/All_isolates.SBS96.tsv"
metadata_file = "mutational_patterns/combined/All_isolates_metadata.tsv"

matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)
metadata = pd.read_csv(metadata_file, sep="\t")

matrix.columns = matrix.columns.str.strip()
metadata.columns = metadata.columns.str.strip()
metadata["sample_id"] = metadata["sample_id"].str.strip()

totals = matrix.sum(axis=0).rename("SBS96_total")
totals.index.name = "sample_id"

# Join metadata to the count totals by sample ID.
metadata = metadata.merge(totals.reset_index(), on="sample_id", how="left", validate="one_to_one")

if metadata["SBS96_total"].isna().any():
    raise ValueError(f"{metadata['SBS96_total'].isna().sum()} samples failed to match SBS96 counts")

print("Overall SBS96 counts per isolate")
print(totals.describe(percentiles=[.05, .25, .50, .75, .95]).to_string())

## since 233 only contains the A and G variants
known = metadata[metadata["MSH6_genotype"].isin(["A233", "G233"])]
print("SBS96 counts by MSH6 genotype")
print(known.groupby("MSH6_genotype")["SBS96_total"].describe(percentiles=[.25, .50, .75]).to_string())
print("Isolates with zero SBS96 counts")
print(int((totals == 0).sum()))
print("Counts by genotype, including unknown")
print(metadata.groupby("MSH6_genotype")["SBS96_total"].agg(["count", "min", "median", "mean", "max"]).to_string())
