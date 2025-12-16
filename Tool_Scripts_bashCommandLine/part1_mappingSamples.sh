#!/usr/bin/env bash
set -euo pipefail

SAMPLEDIR='/projects/ch-lee-lab/USERS/loftum/leonardo/data'
SAMPLEOUTDIRALN='/projects/ch-lee-lab/USERS/loftum/leonardo/aln_mappings'
SAMPLEOUTDIRMEM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_mappings'
REF='/projects/ch-lee-lab/HGSVC/chrY_paper_data/HHU/references/T2T_chm13_122XM_GRCh38-chrY.fasta'

mkdir -p "$SAMPLEOUTDIRALN" "$SAMPLEOUTDIRMEM"

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

OUT_BAMALN="$SAMPLEOUTDIRALN/${FILENAME}.sorted.bam"
OUT_BAMALND="$SAMPLEOUTDIRALN/${FILENAME}.sorted.deduped.bam"
OUT_FLAGSTATALN="$SAMPLEOUTDIRALN/${FILENAME}.flagstat.txt"

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

############################
# 2) BWA-ALN pipeline (legacy)
############################
if [[ ! -f "$OUT_BAMALND" ]]; then
    echo "Running BWA-ALN for $FILENAME..."

    SAI1="$SAMPLEOUTDIRALN/${FILENAME}_R1.sai"
    SAI2="$SAMPLEOUTDIRALN/${FILENAME}_R2.sai"
    TMP_NAME_ALN="$SAMPLEOUTDIRALN/${FILENAME}.namesort.bam"
    TMP_FIX_ALN="$SAMPLEOUTDIRALN/${FILENAME}.fixmate.bam"

    # aln + sampe to unsorted BAM
    bwa aln -t 16 -l 16500 -n 0.01 -o 2 "$REF" "$R1" > "$SAI1"
    bwa aln -t 16 -l 16500 -n 0.01 -o 2 "$REF" "$R2" > "$SAI2"

    bwa sampe "$REF" "$SAI1" "$SAI2" "$R1" "$R2" > "$TMP_NAME_ALN"

    # name-sort (if sampe not already name-sorted; being explicit is safest)
    samtools sort -n -@16 -o "$TMP_NAME_ALN" "$TMP_NAME_ALN"

    # fixmate
    samtools fixmate -m -@16 "$TMP_NAME_ALN" "$TMP_FIX_ALN"

    # coordinate sort
    samtools sort -@16 -o "$OUT_BAMALN" "$TMP_FIX_ALN"

    # mark duplicates
    samtools markdup -@16 "$OUT_BAMALN" "$OUT_BAMALND"

    samtools index "$OUT_BAMALND"
    samtools flagstat "$OUT_BAMALND" > "$OUT_FLAGSTATALN"

    # optional cleanup
    rm -f "$SAI1" "$SAI2" "$TMP_NAME_ALN" "$TMP_FIX_ALN" "$OUT_BAMALN"
else
    echo "BWA-ALN deduped BAM exists for $FILENAME, skipping ALN mapping."
fi
