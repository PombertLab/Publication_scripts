#!/usr/bin/env bash

FASTP=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP

# hellem
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50451,50604,Swiss,50504}.fasta \
  -pe1 $FASTP/{50451,50604,Swiss}.*.R1.fastq.gz \
  -pe2 $FASTP/{50451,50604,Swiss}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# cuniculi
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,GB_M1,ECII_CZ,ECIII_L,ECI,ECII,ECIII}.fasta \
  -pe1 $FASTP/{50602,ECIII_L,ECI,ECII,ECIII}.*.R1.fastq.gz \
  -pe2 $FASTP/{50602,ECIII_L,ECI,ECII,ECIII}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# ECII-CZ => single end mode!
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,GB_M1,ECII_CZ,ECIII_L,ECI,ECII,ECIII}.fasta \
  -fq $FASTP/ECII_CZ.*.SE.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# intestinalis
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50506,50507,50603,50651}.fasta \
  -pe1 $FASTP/{50506,50507,50603,50651}.*.R1.fastq.gz \
  -pe2 $FASTP/{50506,50507,50603,50651}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  --min-avg-qual 30 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA

# Ordospora (no OC4 in SRA)
get_SNPs.pl \
  -fa $GENDIV/FASTA/{OC4,GBEP,NOV7,FISK}.fasta \
  -pe1 $FASTP/{GBEP,NOV7,FISK}.*.R1.fastq.gz \
  -pe2 $FASTP/{GBEP,NOV7,FISK}.*.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  --min-avg-qual 30 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/ILLUMINA
