#!/usr/bin/env bash

#####################################################################
## Downloading data from NCBI Sequence Read Archive (SRA)
## Requires NCBI SRA Toolkit
## https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

DIR=/mnt/d/Linux/PSEU

mkdir -p $DIR
mkdir -p $DIR/CCMP1203
mkdir -p $DIR/UIO007

echo "Downloading dataset Pycnococcus provasolii CCMP1203 from the NCBI SRA archive..."
fasterq-dump SRR28475445 \
  --outdir $DIR/CCMP1203 \
  --progress

echo "Downloading dataset Pseudoscourfieldia marina UIO007 from the NCBI SRA archive..."
fasterq-dump SRR28476106 \
  --outdir $DIR/UIO007 \
  --progress


#####################################################################
## Counting kmers with Jellyfish 2.3.1
## https://github.com/gmarcais/Jellyfish
## On Ubuntu 24.04.3 LTS: sudo apt install jellyfish

## Pycnococcus provasolii CCMP1203
## Pseudoscourfieldia marina UIO007

mkdir -p $DIR/Jellyfish

## Counting kmer + preparing histograms 
for species in {CCMP1203,UIO007}; do

    jellyfish \
    count \
    -C \
    -m 21 \
    -t 16 \
    -s 5G \
    $DIR/$species/*.fastq \
    -o $DIR/Jellyfish/$species.mer

    jellyfish \
    histo \
    -h 3000000 \
    -t 16 \
    -o $DIR/Jellyfish/$species.histo \
    $DIR/Jellyfish/$species.mer

done

#####################################################################
## Plotting histograms (.histo files) manually with GenomeScope 2.0
## https://doi.org/10.1038/s41467-020-14998-3
## http://genomescope.org/genomescope2.0/

