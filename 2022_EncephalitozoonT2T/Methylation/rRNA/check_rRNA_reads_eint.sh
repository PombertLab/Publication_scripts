#!/usr/bin/bash

## Encephalitozoon intestinalis strain ATCC 50506
# Species variables
FASTA=GCA_024399295.1_ASM2439929v1_genomic.fna
PREFIX=eint

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
  "CP075158.1:4223-7994" "CP075158.1:181641-185414" \
  "CP075159.1:4477-8242" "CP075159.1:205462-209228" \
  "CP075160.1:4969-8737" "CP075160.1:208954-212724" \
  "CP075161.1:4750-8516" "CP075161.1:219544-223309" \
  "CP075162.1:4550-8312" "CP075162.1:217954-221725" \
  "CP075163.1:4309-8076" "CP075163.1:222378-226145" \
  "CP075164.1:4532-8301" "CP075164.1:231392-235158" \
  "CP075165.1:4567-8338" "CP075165.1:238771-240054" \
  "CP075166.1:4523-8295" "CP075166.1:260922-264689" \
  "CP075167.1:4398-8169" "CP075167.1:270065-273837" \
  "CP075168.1:4375-8143" "CP075168.1:262093-265856" \
  > $SUB_BAM

samtools bam2fq $SUB_BAM | gzip > $SUB_FASTQ

## Calculating metrics
read_len_plot.py \
  -f $SUB_FASTQ \
  -c deeppink \
  -o ${PREFIX}_rRNA.pdf


