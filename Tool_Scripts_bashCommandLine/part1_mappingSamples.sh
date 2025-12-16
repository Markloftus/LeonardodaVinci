#!/usr/bin/env bash
set -euo pipefail

SAMPLEDIR='/projects/ch-lee-lab/USERS/loftum/leonardo/data'
SAMPLEOUTDIRMEM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_mappings'
REF='/projects/ch-lee-lab/HGSVC/chrY_paper_data/HHU/references/T2T_chm13_122XM_GRCh38-chrY.fasta'

mkdir -p "$SAMPLEOUTDIRMEM"

# Derive sample name from first arg: expects something like SR-0008_human_R1.fastq
in="$1"
FILENAME=$(basename "$in" | awk -F"_human" '{print $1}')

# Original FASTQs
R1_ORIG="$SAMPLEDIR/${FILENAME}_human_R1.fastq"
R2_ORIG="$SAMPLEDIR/${FILENAME}_human_R2.fastq"

# Repaired FASTQs
R1="$SAMPLEDIR/${FILENAME}_human_R1.repair.fastq"
R2="$SAMPLEDIR/${FILENAME}_human_R2.repair.fastq"
SINGLETONS="$SAMPLEDIR/${FILENAME}_human_singletons.fastq"

### REPAIR THE FASTQ FILES (only if not already repaired)
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "Running repair.sh for $FILENAME..."
    repair.sh \
        in1="$R1_ORIG" \
        in2="$R2_ORIG" \
        out1="$R1" \
        out2="$R2" \
        outs="$SINGLETONS"
else
    echo "Repaired FASTQs already exist for $FILENAME, skipping repair."
fi

OUT_BAMMEM="$SAMPLEOUTDIRMEM/${FILENAME}.sorted.bam"
OUT_BAMMEMD="$SAMPLEOUTDIRMEM/${FILENAME}.sorted.deduped.bam"
OUT_FLAGSTATMEM="$SAMPLEOUTDIRMEM/${FILENAME}.flagstat.txt"

############################
# 1) BWA-MEM pipeline
############################
if [[ ! -f "$OUT_BAMMEMD" ]]; then
    echo "Running BWA-MEM for $FILENAME..."

    TMP_NAME="$SAMPLEOUTDIRMEM/${FILENAME}.namesort.bam"
    TMP_FIX="$SAMPLEOUTDIRMEM/${FILENAME}.fixmate.bam"

    # 1) name-sort from bwa mem
    bwa mem -t 16 "$REF" "$R1" "$R2" \
        | samtools sort -n -@16 -o "$TMP_NAME" -

    # 2) add mate tags
    samtools fixmate -m -@16 "$TMP_NAME" "$TMP_FIX"

    # 3) coordinate sort
    samtools sort -@16 -o "$OUT_BAMMEM" "$TMP_FIX"

    # 4) mark duplicates
    samtools markdup -@16 "$OUT_BAMMEM" "$OUT_BAMMEMD"

    # 5) index + stats on deduped BAM
    samtools index "$OUT_BAMMEMD"
    samtools flagstat "$OUT_BAMMEMD" > "$OUT_FLAGSTATMEM"

    # optional cleanup
    rm -f "$TMP_NAME" "$TMP_FIX" "$OUT_BAMMEM"
else
    echo "BWA-MEM deduped BAM exists for $FILENAME, skipping MEM mapping."
fi
