#!/usr/bin/bash

## Encephalitozoon hellem strain ATCC 50604
# Species variables
FASTA=GCA_024399255.1_ASM2439925v1_genomic.fna
PREFIX=ehel

# mapping variables
FASTQ=$PREFIX.nanopore.fastq
MAPPER=winnowmap
CPU=16
BAM="${MAPPER}.BAM/${FASTQ}.${FASTA}.${MAPPER}.bam"
SUB_BAM=${PREFIX}_rRNA.bam
SUB_FASTQ=${PREFIX}_rRNA.fastq.gz

## Mapping
get_SNPs.pl \
  -fa $FASTA \
  -fq $FASTQ \
  -mapper $MAPPER \
  -preset map-ont \
  -t $CPU \
  -rmo \
  -bam

## Creating read subset covering rRNA gene loci with samview
## Must change loci values accordingly; specific to each genome
samtools view \
  -b \
  -h $BAM \
  "CP075147.1:3135-6963" "CP075147.1:252022-255853" \
  "CP075148.1:2666-6495" "CP075148.1:211571-215398" \
  "CP075149.1:3073-6903" "CP075149.1:219144-222973" \
  "CP075150.1:3190-7020" "CP075150.1:237592-241420" \
  "CP075151.1:3164-6989" "CP075151.1:241399-245229" \
  "CP075152.1:3112-6947" "CP075152.1:238433-242261" \
  "CP075153.1:3108-6935" "CP075153.1:246966-250799" \
  "CP075154.1:3537-7365" "CP075154.1:247450-251280" \
  "CP075155.1:3116-6943" "CP075155.1:204463-208293" \
  "CP075156.1:2735-6563" "CP075156.1:270220-274050" \
  "CP075157.1:3063-6891" "CP075157.1:262441-266272" \
  > $SUB_BAM

samtools bam2fq $SUB_BAM | gzip > $SUB_FASTQ

## Calculating metrics
read_len_plot.py \
  -f $SUB_FASTQ \
  -c deeppink \
  -o ${PREFIX}_rRNA.pdf
