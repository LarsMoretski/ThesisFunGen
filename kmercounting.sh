#!/bin/bash

# Author: Lars Moret
#
# Description:
#     Count k-mers for ONE sample using Jellyfish.
#
# Resume behavior:
#     - Valid .jf -> skip
#     - Missing/empty/invalid .jf -> regenerate
#
# Usage:
#     kmercounting.sh --indir /path/to/reads --outdir /path/to/jellyfish --sample SAMPLE_NAME --k 27 --threads 60


set -uo pipefail
THREADS=40
KMER_SIZE=31
HASH_SIZE="2G"


while [[ $# -gt 0 ]]; do
    case "$1" in
        --indir)
            READ_DIR="$2"
            shift 2
            ;;

        --outdir)
            JF_COUNT_DIR="$2"
            shift 2
            ;;

        --sample)
            SAMPLE="$2"
            shift 2
            ;;

        --k)
            KMER_SIZE="$2"
            shift 2
            ;;

        --threads)
            THREADS="$2"
            shift 2
            ;;

        --hash-size)
            HASH_SIZE="$2"
            shift 2
            ;;

        *)
            echo "Unknown argument: $1"
            exit 1
            ;;

    esac

done


if [[ -z "${READ_DIR:-}" ]]; then
    echo "ERROR: --indir is required"
    exit 1
fi

if [[ -z "${JF_COUNT_DIR:-}" ]]; then
    echo "ERROR: --outdir is required"
    exit 1
fi

if [[ -z "${SAMPLE:-}" ]]; then
    echo "ERROR: --sample is required"
    exit 1
fi


mkdir -p "$JF_COUNT_DIR"



valid_jellyfish() {

    local jf_file="$1"

    if [[ ! -f "$jf_file" ]]; then
        return 1
    fi

    if [[ ! -s "$jf_file" ]]; then
        return 1
    fi

    local stats

    if ! stats=$(jellyfish stats "$jf_file" 2>/dev/null); then
        return 1
    fi

    local distinct

    distinct=$(echo "$stats" | awk '/Distinct:/ {print $2; exit}')

    if [[ -z "$distinct" ]]; then
        return 1
    fi

    if [[ "$distinct" -le 0 ]]; then
        return 1
    fi

    return 0
}


R1="${READ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${READ_DIR}/${SAMPLE}_R2.fastq.gz"

if [[ ! -f "$R1" ]]; then
    echo "$R1 missing R1"
    exit 1

fi

if [[ ! -f "$R2" ]]; then
    echo "$R2 missing R2"
    exit 1

fi


OUTFILE="${JF_COUNT_DIR}/${SAMPLE}.jf"
TMPFILE="${OUTFILE}.tmp.$$"

if valid_jellyfish "$OUTFILE"; then
    echo "${SAMPLE}.jf is valid"
    exit 0

fi

if [[ -e "$OUTFILE" ]]; then
    echo "Existing .jf is invalid"
    rm -f "$OUTFILE"

else
    echo "Jellyfish database does not exist"
fi

rm -f "$TMPFILE"

if jellyfish count \
    -m "$KMER_SIZE" \
    -C \
    -s "$HASH_SIZE" \
    -t "$THREADS" \
    <(zcat "$R1" "$R2") \
    -o "$TMPFILE"
then

    if valid_jellyfish "$TMPFILE"; then
        mv "$TMPFILE" "$OUTFILE"
        echo "Finished ${SAMPLE}"
        exit 0

    else
        echo "Jellyfish produced an invalid database"
        rm -f "$TMPFILE"
        exit 1
    fi

else
    echo "Jellyfish failed for ${SAMPLE}"
    rm -f "$TMPFILE"
    exit 1

fi
