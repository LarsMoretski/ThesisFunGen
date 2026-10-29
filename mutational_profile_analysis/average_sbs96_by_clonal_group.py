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

## import statements
from pathlib import Path
import sys
import numpy as np
import pandas as pd


BASE = Path("/lustre/BIF/nobackup/moret013/thesis_GEN")
INPUT_MATRIX = (BASE / "mutational_patterns/combined/All_isolates.SBS96_proportions.tsv")
INPUT_METADATA = (BASE / "mutational_patterns/combined/All_isolates_metadata.tsv")
OUTPUT_DIR = BASE / "mutational_patterns/combined"
OUTPUT_MATRIX = OUTPUT_DIR / "ClonalGroup_SBS96_mean_proportions.tsv"
OUTPUT_METADATA = OUTPUT_DIR / "ClonalGroup_SBS96_group_metadata.tsv"
EXCLUDED_GROUPS = {"ClonalGroup_002"} ##Clonal group 2 was excluded due to haplotype mixture
VALID_GENOTYPES = {"A233", "G233"}  ## There are no genotypes recorded other than A233 and G233


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(INPUT_MATRIX, sep="\t", index_col=0)
    metadata = pd.read_csv(INPUT_METADATA, sep="\t", dtype=str)
    required_columns = {"sample_id", "clonal_group", "MSH6_genotype"}
    missing = required_columns - set(metadata.columns)
    if missing:
        sys.exit(f"ERROR: Metadata is missing columns: {sorted(missing)}")

    if matrix.index.has_duplicates:
        sys.exit("ERROR: Duplicate SBS96 context labels in matrix.")

    if matrix.columns.has_duplicates:
        sys.exit("ERROR: Duplicate sample IDs in SBS96 matrix.")

    if metadata["sample_id"].duplicated().any():
        sys.exit("ERROR: Duplicate sample IDs in metadata.")

    matrix = matrix.apply(pd.to_numeric, errors="raise")
    values = matrix.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        sys.exit("ERROR: Matrix contains missing or non-finite values.")

    if (values < 0).any():
        sys.exit("ERROR: Matrix contains negative proportions.")

    if not np.allclose(values.sum(axis=0), 1.0, atol=1e-6):
        sys.exit("ERROR: Some isolate spectra do not sum to 1.")

    metadata = metadata.set_index("sample_id", drop=False)
    matrix_samples = set(matrix.columns)
    metadata_samples = set(metadata.index)
    missing_metadata = matrix_samples - metadata_samples
    missing_matrix = metadata_samples - matrix_samples

    if missing_metadata:
        sys.exit(
            f"ERROR: {len(missing_metadata)} matrix samples lack metadata. ")

    if missing_matrix:
        print(
            f"NOTE: {len(missing_matrix)} metadata samples are not in the SBS96 matrix and won't be used")

    metadata = metadata.loc[list(matrix.columns)].copy()   # Keep only A233/G233 isolates, excluding the mixed clonal group.
    group_status = metadata.groupby("clonal_group").agg(n_isolates=("sample_id", "size"), n_genotypes=("MSH6_genotype", "nunique"),
        all_known=("MSH6_genotype", lambda x: x.isin(VALID_GENOTYPES).all()),)

    eligible_groups = group_status.index[(group_status["n_genotypes"] == 1) & group_status["all_known"] & ~group_status.index.isin(EXCLUDED_GROUPS)]
    selected_metadata = metadata[metadata["clonal_group"].isin(eligible_groups)].copy()

    if selected_metadata.empty:
        sys.exit("ERROR: No eligible isolates remain after filtering.")

    genotype_counts = selected_metadata.groupby("clonal_group")["MSH6_genotype"].nunique()
    mixed_groups = genotype_counts[genotype_counts > 1]
    if not mixed_groups.empty:
        sys.exit(
            "ERROR: Mixed-genotype groups remain after filtering: "
            + ", ".join(mixed_groups.index))

    # One mean spectrum per clonal group so that every isolate contributes equally
    selected_samples = selected_metadata.index.tolist()
    selected_matrix = matrix.loc[:, selected_samples]

    group_profiles = []
    group_records = []
    for group, group_meta in selected_metadata.groupby("clonal_group", sort=True):
        sample_ids = group_meta.index.tolist()
        genotype = group_meta["MSH6_genotype"].iloc[0]

        # Average each of the 96 context proportions across isolates.
        mean_profile = selected_matrix.loc[:, sample_ids].mean(axis=1)
        if not np.isclose(mean_profile.sum(), 1.0, atol=1e-6):
            sys.exit(f"ERROR: Mean spectrum for {group} does not sum to 1 ")

        group_profiles.append(pd.Series(mean_profile.to_numpy(), index=matrix.index, name=group))

        group_records.append({
                "clonal_group": group,
                "MSH6_genotype": genotype,
                "group_size": len(sample_ids),
                "sample_ids": ";".join(sample_ids),})

    profiles_df = pd.DataFrame(group_profiles)
    profiles_df.index.name = "clonal_group"
    profiles_df.to_csv(OUTPUT_MATRIX, sep="\t", float_format="%.10g")

    group_metadata = pd.DataFrame(group_records)
    group_metadata = group_metadata.sort_values(["MSH6_genotype", "clonal_group"])
    group_metadata.to_csv(OUTPUT_METADATA, sep="\t", index=False)

    print(f"Input isolates in matrix: {matrix.shape[1]}")
    print(f"Eligible isolates retained: {len(selected_samples)}")
    print(f"Excluded mixed group(s): {', '.join(sorted(EXCLUDED_GROUPS))}")
    print(f"Groups retained: {len(group_metadata)}")
    print(group_metadata["MSH6_genotype"].value_counts().to_string())


## Main
if __name__ == "__main__":
    main()

