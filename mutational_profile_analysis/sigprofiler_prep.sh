#!/bin/bash
##Author: Lars Moret, MSc student @ FunGen
##Date: 29/10/2026
##
##Description:
##        prepare the VCFs for sigprofiler
##Usage: bash sigprofiler_prep.sh

mkdir -p scripts logs
cat > scripts/split_all_clonal_groups.sh <<'BASH'
set -Eeuo pipefail

##set paths
BASE_DIR="/lustre/BIF/nobackup/moret013/thesis_GEN"
INPUT_DIR="${BASE_DIR}/snelders_pipeline/missing_group_vcfs"
OUTPUT_DIR="${BASE_DIR}/mutational_patterns/all_group_vcfs"
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

##Create logfile
LOG_FILE="${LOG_DIR}/split_all_clonal_groups_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Check dependencies
for cmd in bcftools; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "$cmd is not available in PATH."
        exit 1
    fi
done

shopt -s nullglob
GROUP_FILES=("$INPUT_DIR"/ClonalGroup_*.vcf.gz)

if (( ${#GROUP_FILES[@]} == 0 )); then
    echo "No clonal-group VCFs found in $INPUT_DIR"
    exit 1
fi

echo "Found ${#GROUP_FILES[@]} clonal-group VCFs."

GROUPS_DONE=0
SAMPLES_DONE=0
ERRORS=0

for GROUP_VCF in "${GROUP_FILES[@]}"; do
    GROUP=$(basename "$GROUP_VCF" .vcf.gz)
    GROUP_DIR="${OUTPUT_DIR}/${GROUP}"
    mkdir -p "$GROUP_DIR"
    echo "[$((GROUPS_DONE + 1))/${#GROUP_FILES[@]}] Processing $GROUP"
    SAMPLE_LIST=$(mktemp)
    if ! bcftools query -l "$GROUP_VCF" > "$SAMPLE_LIST"; then
        echo "Cant read sample names from $GROUP_VCF"
        rm -f "$SAMPLE_LIST"
        ERRORS=$((ERRORS + 1))
        continue
    fi

    GROUP_SAMPLES=0

    while IFS= read -r SAMPLE; do
        [[ -z "$SAMPLE" ]] && continue

        OUT_VCF="${GROUP_DIR}/${SAMPLE}.vcf.gz"

        # Skip an existing, non-empty output with a readable index.
        if [[ -s "$OUT_VCF" && -s "${OUT_VCF}.csi" ]]; then
            if bcftools index -n "$OUT_VCF" >/dev/null 2>&1; then
                GROUP_SAMPLES=$((GROUP_SAMPLES + 1))
                continue
            fi
        fi

        TMP_VCF="${OUT_VCF%.vcf.gz}.tmp.$$.vcf.gz"

        if bcftools view \
            -s "$SAMPLE" \
            -c 1 \
            -Oz \
            -o "$TMP_VCF" \
            "$GROUP_VCF" \
        && bcftools index -f "$TMP_VCF"; then

            mv "$TMP_VCF" "$OUT_VCF"
            mv "${TMP_VCF}.csi" "${OUT_VCF}.csi"
            GROUP_SAMPLES=$((GROUP_SAMPLES + 1))
        else
            echo "Failed to extract $SAMPLE from $GROUP"
            rm -f "$TMP_VCF" "${TMP_VCF}.csi"
            ERRORS=$((ERRORS + 1))
        fi

    done < "$SAMPLE_LIST"
    rm -f "$SAMPLE_LIST"

    GROUPS_DONE=$((GROUPS_DONE + 1))
    SAMPLES_DONE=$((SAMPLES_DONE + GROUP_SAMPLES))
done

if (( ERRORS > 0 )); then
    echo "Some files failed, check log before continueing"
    exit 1
fi

echo "All groups processed successfully."




