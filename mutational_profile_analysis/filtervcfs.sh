#!/bin/bash
##Author: Lars Moret, MSc student @ FunGen
##Date: 29/10/2026

#Description:
##    BASH script for filtering VCF files using bcftools to only retain genotypes that differ between samples (clonal groups)
##    Missing genomes are excluded from the comparison

## Usage: bash filtervcfs.sh


set -euo pipefail

INPUT_DIR="snelders_pipeline"
OUTPUT_DIR="clonal_group_variable_vcfs"

mkdir -p "$OUTPUT_DIR"

for vcf in "$INPUT_DIR"/ClonalGroup_*.vcf.gz; do
    [[ -e "$vcf" ]] || continue

    name=$(basename "$vcf" .vcf.gz)
    out="$OUTPUT_DIR/${name}.variable.vcf.gz"

    bcftools view -i 'COUNT(GT="mis")<N_SAMPLES && N_ALT>=0' "$vcf" -Ou \
      | bcftools +fill-tags -Ou -- -t AC,AN \
      | bcftools view -i 'AC>0 && AC<AN' -Oz -o "$out"

    bcftools index -t "$out"
done


