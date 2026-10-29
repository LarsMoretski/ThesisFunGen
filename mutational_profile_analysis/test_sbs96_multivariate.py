#!/usr/bin/env python3
'''
Author: Lars Moret, MSc student @ FunGen
Date: 29/10/2026

Description:
    Python script for an exploratory multivariate analysis

Usage: python3 test_sbs96_multivariate.py
'''

## import packages
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway

## set path
BASE = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
MATRIX = BASE / "mutational_patterns/combined/ClonalGroup_SBS96_mean_proportions.tsv"
META = BASE / "mutational_patterns/combined/ClonalGroup_SBS96_group_metadata.tsv"
OUT = BASE / "mutational_patterns/combined"

## statistical test metrics
N_PERM = 9999
SEED = 42
rng = np.random.default_rng(SEED)

## load files
profiles = pd.read_csv(MATRIX, sep="\t")
meta = pd.read_csv(META, sep="\t")
required = {"clonal_group", "MSH6_genotype", "group_size"}
if not required.issubset(meta.columns):
    raise ValueError(f"Metadata missing columns, {required - set(meta.columns)}")

# Identify the group-ID column in the profile matrix.
group_ids = set(meta["clonal_group"].astype(str))
overlap = {col: profiles[col].astype(str).isin(group_ids).sum() for col in profiles.columns}
id_col = max(overlap, key=overlap.get)
if overlap[id_col] == 0:
    raise ValueError("Could not identify group IDs in the profile matrix")

profiles = profiles.rename(columns={id_col: "clonal_group"})
profiles["clonal_group"] = profiles["clonal_group"].astype(str)
meta["clonal_group"] = meta["clonal_group"].astype(str)
meta = meta[meta["MSH6_genotype"].isin(["A233", "G233"])].copy() ## only keep groups that only contain G233 or A233
if meta["clonal_group"].duplicated().any():
    raise ValueError("Duplicate clonal groups in metadata.")

data = meta.merge(profiles, on="clonal_group", how="inner", validate="one_to_one")

context_cols = [c for c in profiles.columns if c != "clonal_group" and pd.api.types.is_numeric_dtype(profiles[c])]

if len(context_cols) != 96:
    raise ValueError(f"Expected 96 numeric SBS contexts but found {len(context_cols)}")

data[context_cols] = data[context_cols].apply(pd.to_numeric, errors="raise")
X = data[context_cols].to_numpy(dtype=float)
groups = data["MSH6_genotype"].to_numpy()

if not np.isfinite(X).all() or (X < 0).any():
    raise ValueError("Profiles contain missing, infinite, or negative values.")

row_sums = X.sum(axis=1)
if not np.allclose(row_sums, 1, atol=1e-5):
    raise ValueError("Each group profile should sum to 1.")

counts = pd.Series(groups).value_counts()
if not {"A233", "G233"}.issubset(counts.index):
    raise ValueError("Both genotype groups must be represent the same genotype")


Y = np.sqrt(X) ## helling transformation
D = squareform(pdist(Y, metric="euclidean")) ## compute hellinger distance
n = len(groups)

# PERMANOVA pseudo-F from the distance matrix.
def permanova_stat(D, labels):
    n = len(labels)
    unique = np.unique(labels)
    if len(unique) != 2:
        raise ValueError("This implementation expects two groups.")

    d2 = D ** 2
    ss_total = np.triu(d2, 1).sum() / n
    group_ss = 0.0

    for label in unique:
        idx = np.flatnonzero(labels == label)
        group_ss += np.triu(d2[np.ix_(idx, idx)], 1).sum() / len(idx)

    ss_within = group_ss
    ss_between = ss_total - ss_within
    df_between = len(unique) - 1
    df_within = n - len(unique)

    pseudo_f = (ss_between / df_between) / (ss_within / df_within)
    r2 = ss_between / ss_total
    return pseudo_f, r2

observed_f, observed_r2 = permanova_stat(D, groups)
perm_f = np.empty(N_PERM)

for i in range(N_PERM):
    perm_f[i], _ = permanova_stat(D, rng.permutation(groups))

p_permanova = (1 + np.sum(perm_f >= observed_f)) / (N_PERM + 1)
J = np.eye(n) - np.ones((n, n)) / n
B = -0.5 * J @ (D ** 2) @ J
eigvals, eigvecs = np.linalg.eigh(B)
keep = eigvals > 1e-10 ## negative eigenvalues are discarded

if not keep.any():
    raise ValueError("No positive PCoA eigenvalues.")

coords = eigvecs[:, keep] * np.sqrt(eigvals[keep])

def centroid_distances(coords, labels):
    distances = np.empty(len(labels))
    for label in np.unique(labels):
        idx = np.flatnonzero(labels == label)
        centroid = coords[idx].mean(axis=0)
        distances[idx] = np.linalg.norm(coords[idx] - centroid, axis=1)

    return distances

def dispersion_f(distances, labels):
    vals = [distances[labels == g] for g in np.unique(labels)]

    return f_oneway(*vals).statistic

dist_centroid = centroid_distances(coords, groups)
observed_disp_f = dispersion_f(dist_centroid, groups)
perm_disp_f = np.empty(N_PERM)

for i in range(N_PERM):
    permuted = rng.permutation(groups)
    perm_disp_f[i] = dispersion_f(dist_centroid, permuted)

p_disp = (1 + np.sum(perm_disp_f >= observed_disp_f)) / (N_PERM + 1)

results = pd.DataFrame([{
        "test": "PERMANOVA_Hellinger_Euclidean",
        "statistic": observed_f,
        "R2": observed_r2,
        "permutations": N_PERM,
        "p_value": p_permanova,
        "permutation_scheme": "unrestricted",},
    {
        "test": "PERMDISP_PCoA_centroid",
        "statistic": observed_disp_f,
        "R2": np.nan,
        "permutations": N_PERM,
        "p_value": p_disp,
        "permutation_scheme": "unrestricted",}])

results.to_csv(OUT / "SBS96_multivariate_test.tsv", sep="\t", index=False)

pd.DataFrame({
    "clonal_group": data["clonal_group"],
    "MSH6_genotype": groups,
    "distance_to_centroid": dist_centroid}).to_csv(OUT / "SBS96_group_dispersion.tsv", sep="\t", index=False)

