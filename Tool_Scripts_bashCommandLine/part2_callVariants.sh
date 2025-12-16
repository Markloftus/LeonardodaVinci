#!/usr/bin/env bash
set -euo pipefail

POZNIK_COORD='/projects/ch-lee-lab/USERS/loftum/leonardo/scripts/part2_variant/Poznik_2013_S1b.chrY.callable.hg38.bed'
SAMPLEDIR='/projects/ch-lee-lab/USERS/loftum/leonardo/data'
SAMPLEOUTDIRMEMBAM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_mappings'
SAMPLEOUTDIRMEM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_variants'
REF='/projects/ch-lee-lab/HGSVC/chrY_paper_data/HHU/references/T2T_chm13_122XM_GRCh38-chrY.fasta'

in="$1"
FILENAME=$(basename "$in" | awk -F".sorted." '{print $1}')
FILE=$(basename "$in")


#SAMTOOLS
#samtools coverage -o "$SAMPLEOUTDIRMEMBAM/$FILENAME.coverage.txt" "$SAMPLEOUTDIRMEMBAM/$FILE"
#samtools depth -a "$SAMPLEOUTDIRMEMBAM/$FILE" > "$SAMPLEOUTDIRMEMBAM/$FILENAME.depth.txt"

bcftools mpileup \
  --threads 16 \
  -R "$POZNIK_COORD" \
  -f "$REF" \
  -C 0 \
  -q 5 \
  -Q 20 \
  -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR \
  "$SAMPLEOUTDIRMEMBAM/$FILE" \
  -Ou \
| bcftools call \
    --threads 16 \
    -m \
    -Oz \
    -o "$SAMPLEOUTDIRMEM/$FILENAME.map5.vcf.gz"


bcftools index -t "$SAMPLEOUTDIRMEM/$FILENAME.map5.vcf.gz"

bcftools filter \
  -g 5 \
  -Ov \
  -o "$SAMPLEOUTDIRMEM/$FILENAME.map5.sg5.vcf.gz" \
  "$SAMPLEOUTDIRMEM/$FILENAME.map5.vcf.gz"

#indels get removed with python after
#bcftools view \
#  -Ov \
#  -o "$SAMPLEOUTDIRMEM/$FILENAME.sg5.woINDEL.vcf.gz" \
#  "$SAMPLEOUTDIRMEM/$FILENAME.sg5.vcf.gz"
