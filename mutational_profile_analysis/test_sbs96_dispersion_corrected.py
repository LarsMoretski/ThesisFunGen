#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script for computing the pairwise distance between two k-mer indexes using:
        - Jaccard distance
        - Weighted Jaccard similarity
        - Dice-Sørensen coefficient

Usage: python3 jellyfish_to_parquet.py --indir /path/to/files.jf --outdir /path/to/outdir
'''

## import packages
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import f_oneway, mannwhitneyu

## set paths
BASE = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
MATRIX = BASE / "mutational_patterns/combined/ClonalGroup_SBS96_mean_proportions.tsv"
META = BASE / "mutational_patterns/combined/ClonalGroup_SBS96_group_metadata.tsv"
OUT = BASE / "mutational_patterns/combined"

## set metrics
N_PERM = 9999
SEED = 42
rng = np.random.default_rng(SEED)

## Load grouplevel SBS96 profiles and metadata
profiles = pd.read_csv(MATRIX, sep="\t")
meta = pd.read_csv(META, sep="\t")

group_ids = set(meta["clonal_group"].astype(str))
matches = { col: profiles[col].astype(str).isin(group_ids).sum() for col in profiles.columns}
id_col = max(matches, key=matches.get)
if matches[id_col] == 0:
    raise ValueError("Can't identify clonal group IDs in the profile matrix.")

profiles = profiles.rename(columns={id_col: "clonal_group"})
profiles["clonal_group"] = profiles["clonal_group"].astype(str)
meta["clonal_group"] = meta["clonal_group"].astype(str)
meta = meta[meta["MSH6_genotype"].isin(["A233", "G233"])].copy()

if meta["clonal_group"].duplicated().any():
    raise ValueError("Duplicate group IDs in metadata.")

data = meta.merge(profiles, on="clonal_group", how="inner", validate="one_to_one")
context_cols = [c for c in profiles.columns if c != "clonal_group" and pd.api.types.is_numeric_dtype(profiles[c])]

if len(context_cols) != 96:
    raise ValueError(f"Expected 96 SBS files, found {len(context_cols)}")

X = data[context_cols].to_numpy(dtype=float)
labels = data["MSH6_genotype"].to_numpy()

if not np.isfinite(X).all() or (X < 0).any():
    raise ValueError("Profiles contain missing, infinite, or negative values.")

if not np.allclose(X.sum(axis=1), 1, atol=1e-5):
    raise ValueError("At least one SBS96 profile does not sum to 1.")

if set(labels) != {"A233", "G233"}:
    raise ValueError("Both A233 and G233 must be present.")

Y = np.sqrt(X) ## hellinger transformation

def distances_to_centroids(coords, groups):
    """Recalculate each group's centroid and distances for these labels."""
    distances = np.empty(len(groups), dtype=float)
    for group in np.unique(groups):
        idx = np.flatnonzero(groups == group)
        centroid = coords[idx].mean(axis=0)
        distances[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)

    return distances



def dispersion_f(distances, groups):
    values = [distances[groups == group] for group in np.unique(groups)]

    return float(f_oneway(*values).statistic)


observed_distances = distances_to_centroids(Y, labels)
observed_f = dispersion_f(observed_distances, labels)

perm_f = np.empty(N_PERM) ##recompute centroids and distances to the centroids for every permutation
for i in range(N_PERM):
    permuted_labels = rng.permutation(labels)
    permuted_distances = distances_to_centroids(Y, permuted_labels)
    perm_f[i] = dispersion_f(permuted_distances, permuted_labels)

p_value = (1 + np.sum(perm_f >= observed_f)) / (N_PERM + 1)

summary = (pd.DataFrame({"MSH6_genotype": labels, "distance_to_centroid": observed_distances}).groupby("MSH6_genotype")["distance_to_centroid"].agg(["count", "mean", "std", "median", "min", "max"]))

mw = mannwhitneyu(observed_distances[labels == "A233"], observed_distances[labels == "G233"], alternative="two-sided")

result = pd.DataFrame([{
    "test": "Corrected_permutation_dispersion_test",
    "distance": "Euclidean_on_Hellinger_profiles",
    "statistic_F": observed_f,
    "permutations": N_PERM,
    "p_value": p_value,
    "seed": SEED,
    "permutation_scheme": "unrestricted_labels_centroids_recomputed"}])

result.to_csv(OUT / "SBS96_dispersion_corrected_test.tsv", sep="\t", index=False)
summary.to_csv(OUT / "SBS96_dispersion_corrected_summary.tsv", sep="\t")

pd.DataFrame({"clonal_group": data["clonal_group"], "MSH6_genotype": labels, "distance_to_centroid": observed_distances}).to_csv(OUT / "SBS96_dispersion_corrected_distances.tsv", sep="\t", index=False)
print(f"{len(labels)} groups analyzed")
print("Observed distances to genotype centroid:")
print(summary.to_string())
print("Descriptive Mann-Whitney test:")
print(f"U = {mw.statistic:.3f}, p = {mw.pvalue:.6g}")

print("Corrected permutation dispersion test:")
print(result.to_string(index=False))
print("Saved corrected results under", OUT)
