#!/bin/bash
#"""
#   Author: Lars Moret
#   Description: Bash script for running the K-mer counting analysis using jellyfish
#                Value for K is optimized based on a dataset of ~2k A. fumigatus isolates
#   Usage:
#"""

set -euo pipefail
## reports and aborts when "silent" errors occur, instead of ignoring

## The directories are hardcoded in this script. Changing this will corrupt
## the NextFlow pipeline, so do at your own risk.

read_dir="/lustre/BIF/nobackup/moret013/thesis_GEN/final_subdata/subsamples/10000_reads"
sample_list="/lustre/BIF/nobackup/moret013/thesis_GEN/subdataset/filtered_accession.txt"
jf_count_dir="/lustre/BIF/nobackup/moret013/thesis_GEN/final_subdata/jellyfish/10k_31/"
jf_dump_dir="/lustre/BIF/nobackup/moret013/thesis_GEN/subdataset/jellyfish/10k_31_dump"

## alt_kmer_size is discontinued after testing
alt_kmer_size=51
## change the kmer size, hash size and n threads here (NOT IN THE SCRIPT!)
kmer_size=31
hash_size="2G"
threads=40

## create the output dir if it doesnt exist
mkdir -p "$jf_count_dir"
mkdir -p "$jf_dump_dir"

while read -r sample
do
    echo "checking if ${sample} exists"
    R1="${read_dir}/${sample}_R1.fastq.gz"
    R2="${read_dir}/${sample}_R2.fastq.gz"

    if [[ ! -f "$R1" ]]; then
        echo "seems like $R1 is missing"
        continue
    fi
    if [[ ! -f "$R2" ]]; then
        echo "seems like $R2 is missing"
        continue
    fi

    jellyfish count -m ${kmer_size} -C -s ${hash_size} -t ${threads} <(zcat "$R1" "$R2") -o "${jf_count_dir}/${sample}.jf"
    echo "finished kmer counting ${sample}"

done < "$sample_list"
echo "kmer counting of all samples in ${sample_list} finished, starting kmer extraction"

for jf in ${jf_count_dir}/*.jf
do
    sample2=$(basename "$jf" .jf)
    echo "dumping ${sample2}"
    jellyfish dump -c "$jf" | sort | gzip > "${jf_dump_dir}/${sample2}.kmercounts.gz"

done
echo "finished extracting kmers"
