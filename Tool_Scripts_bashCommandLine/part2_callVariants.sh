#!/usr/bin/env bash
set -euo pipefail

POZNIK_COORD='/projects/ch-lee-lab/USERS/loftum/leonardo/scripts/part2_variant/Poznik_2013_S1b.chrY.callable.hg38.bed'
SAMPLEDIR='/projects/ch-lee-lab/USERS/loftum/leonardo/data'
SAMPLEOUTDIRALNBAM='/projects/ch-lee-lab/USERS/loftum/leonardo/aln_mappings'
SAMPLEOUTDIRMEMBAM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_mappings'
SAMPLEOUTDIRALN='/projects/ch-lee-lab/USERS/loftum/leonardo/aln_variants'
SAMPLEOUTDIRMEM='/projects/ch-lee-lab/USERS/loftum/leonardo/mem_variants'
REF='/projects/ch-lee-lab/HGSVC/chrY_paper_data/HHU/references/T2T_chm13_122XM_GRCh38-chrY.fasta'

in="$1"
FILENAME=$(basename "$in" | awk -F".sorted." '{print $1}')
FILE=$(basename "$in")


#SAMTOOLS
#samtools coverage -o "$SAMPLEOUTDIRALNBAM/$FILENAME.coverage.txt" "$SAMPLEOUTDIRALNBAM/$FILE"
#samtools coverage -o "$SAMPLEOUTDIRMEMBAM/$FILENAME.coverage.txt" "$SAMPLEOUTDIRMEMBAM/$FILE"

#samtools depth -a "$SAMPLEOUTDIRALNBAM/$FILE" > "$SAMPLEOUTDIRALNBAM/$FILENAME.depth.txt"
#samtools depth -a "$SAMPLEOUTDIRMEMBAM/$FILE" > "$SAMPLEOUTDIRMEMBAM/$FILENAME.depth.txt"

bcftools mpileup \
  --threads 16 \
  -R "$POZNIK_COORD" \
  -f "$REF" \
  -C 0 \
  -q 5 \
  -Q 20 \
  -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD,INFO/ADF,INFO/ADR \
  "$SAMPLEOUTDIRALNBAM/$FILE" \
  -Ou \
| bcftools call \
    --threads 16 \
    -m \
    -Oz \
    -o "$SAMPLEOUTDIRALN/$FILENAME.map5.vcf.gz"


bcftools index -t "$SAMPLEOUTDIRALN/$FILENAME.map5.vcf.gz"

bcftools filter \
  -g 5 \
  -Ov \
  -o "$SAMPLEOUTDIRALN/$FILENAME.map5.sg5.vcf.gz" \
  "$SAMPLEOUTDIRALN/$FILENAME.map5.vcf.gz"

#bcftools view \
#  -v snps \
#  -Ov \
#  -o "$SAMPLEOUTDIRALN/$FILENAME.sg5.woINDEL.vcf.gz" \
#  "$SAMPLEOUTDIRALN/$FILENAME.sg5.vcf.gz"

##############################################
#DO THE SAME FOR THE MEM MAPPINGS
##############################################

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

#bcftools view \
#  -v snps \
#  -Ov \
#  -o "$SAMPLEOUTDIRMEM/$FILENAME.sg5.woINDEL.vcf.gz" \
#  "$SAMPLEOUTDIRMEM/$FILENAME.sg5.vcf.gz"
