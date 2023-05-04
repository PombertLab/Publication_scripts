#!/usr/bin/bash

## Encephalitozoon cuniculi strain ATCC 50602
# Species variables
FASTA=ecun.fasta
PREFIX=ecun

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
  "CP091431:5362-12772" "CP091431:290602-298013" \
  "CP091432:5209-12620" "CP091432:287148-294558" \
  "CP091433:5436-12848" "CP091433:274382-281792" \
  "CP091434:5320-12729" "CP091434:257143-264554" \
  "CP091435:4966-12377" "CP091435:248680-256090" \
  "CP091436:5288-12699" "CP091436:236021-243432" \
  "CP091437:5031-12441" "CP091437:233431-240842" \
  "CP091438:5065-12475" "CP091438:230117-237530" \
  "CP091439:5445-12856" "CP091439:228046-235456" \
  "CP091440:4948-12358" "CP091440:221259-228719" \
  "CP091441:5298-12714" "CP091441:203953-211362" \
  > $SUB_BAM

samtools bam2fq $SUB_BAM | gzip > $SUB_FASTQ

## Calculating metrics
read_len_plot.py \
  -f $SUB_FASTQ \
  -c deeppink \
  -o ${PREFIX}_rRNA.pdf


