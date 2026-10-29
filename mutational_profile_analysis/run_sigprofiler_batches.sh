#!/bin/bash
##Author: Lars Moret, MSc student @ FunGen
##Date: 29/10/2026
##
##Description:
##    Python script for computing the pairwise distance between two k-mer indexes using:
        - Jaccard distance
        - Weighted Jaccard similarity
        - Dice-Sørensen coefficient

##Usage: bash run_sigprofiler_batches.sh

set -uo pipefail

## set path
ROOT="/lustre/BIF/nobackup/moret013/thesis_GEN"
SOURCE="${ROOT}/mutational_patterns/all_group_vcfs"
WORK="${ROOT}/mutational_patterns/full_batches"

BATCH_SIZE=50

## check environment
if ! python -c \
  "from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc" \
  >/dev/null 2>&1; then
    echo "SigProfiler is not available in this Python environment."
    exit 1
fi

mkdir -p "$WORK" "$ROOT/logs"

MANIFEST="${WORK}/isolate_vcfs.txt"
if [[ ! -s "$MANIFEST" ]]; then
    find "$SOURCE" -type f -name '*.vcf.gz' | sort > "$MANIFEST"
fi

TOTAL=$(wc -l < "$MANIFEST")
BATCHES=$(( (TOTAL + BATCH_SIZE - 1) / BATCH_SIZE ))

if [[ "$TOTAL" -eq 0 ]]; then
    echo "No VCFs found in $SOURCE"
    exit 1
fi

FAILED=0

for ((BATCH=1; BATCH<=BATCHES; BATCH++)); do
    printf -v PROJECT "Batch%03d" "$BATCH"
    BATCH_DIR="${WORK}/${PROJECT}"
    INPUT_DIR="${BATCH_DIR}/input"
    OUTPUT_DIR="${BATCH_DIR}/results"
    LOG="${BATCH_DIR}/${PROJECT}.log"
    DONE="${BATCH_DIR}/COMPLETE"

    START=$(( (BATCH - 1) * BATCH_SIZE + 1 ))
    END=$(( BATCH * BATCH_SIZE ))
    if (( END > TOTAL )); then END=$TOTAL; fi

    mkdir -p "$BATCH_DIR"

    if [[ -f "$DONE" &&
          -s "${OUTPUT_DIR}/SBS/${PROJECT}.SBS96.all" ]]; then
        echo "[$BATCH/$BATCHES] $PROJECT already complete; skipping."
        continue
    fi

    echo "[$BATCH/$BATCHES] Processing $PROJECT (isolates $START-$END)"

    rm -rf "$INPUT_DIR"
    mkdir -p "$INPUT_DIR"
    rm -f "$DONE"

    BATCH_LIST="${BATCH_DIR}/vcfs.txt"
    sed -n "${START},${END}p" "$MANIFEST" > "$BATCH_LIST"

    OK=1
    while IFS= read -r VCF; do
        SAMPLE=$(basename "$VCF" .vcf.gz)
        if ! gzip -cd "$VCF" > "${INPUT_DIR}/${SAMPLE}.vcf"; then
            echo "Failed to decompress $VCF"
            OK=0
            break
        fi
    done < "$BATCH_LIST"

    if [[ "$OK" -ne 1 ]]; then
        echo "FAILED: $PROJECT (input prep)"
        FAILED=$((FAILED + 1))
        continue
    fi

    EXPECTED=$(wc -l < "$BATCH_LIST")
    PREPARED=$(find "$INPUT_DIR" -maxdepth 1 -type f \
        -name '*.vcf' | wc -l)

    if [[ "$PREPARED" -ne "$EXPECTED" ]]; then
        echo "$PROJECT failed, expected $EXPECTED inputs but found $PREPARED"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Run SigProfiler in python and keep a separate log for each batch
    ### !!!! note that the reference genome has to be fitted manually and has to pass MSD5 checksum
    ### therefore, that name can change, so check if your reference genome is also named Af293!!!
    if PROJECT="$PROJECT" INPUT_DIR="$INPUT_DIR" \
       OUTPUT_DIR="$OUTPUT_DIR" \
       python - <<'PY' > "$LOG" 2>&1

import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc

SigProfilerMatrixGeneratorFunc.SigProfilerMatrixGeneratorFunc(
    project=os.environ["PROJECT"],
    reference_genome="Af293",
    path_to_input_files=os.environ["INPUT_DIR"] + "/",
    output_directory=os.environ["OUTPUT_DIR"] + "/",
    exome=False,
    chrom_based=False,
    plot=True,
    tsb_stat=True,
    seqInfo=True
)
PY
    then
        if [[ -s "${OUTPUT_DIR}/SBS/${PROJECT}.SBS96.all" ]]; then
            date > "$DONE"
            echo "$PROJECT passed ($EXPECTED isolates)"
        else
            echo "$PROJECT has failed. SBS96 matrix missing, see $LOG for more information"
            FAILED=$((FAILED + 1))
        fi
    else
        echo "$PROJECT has failed(SigProfiler error, see $LOG for more information)"
        FAILED=$((FAILED + 1))
    fi
done

if [[ "$FAILED" -gt 0 ]]; then
    exit 1
fi
