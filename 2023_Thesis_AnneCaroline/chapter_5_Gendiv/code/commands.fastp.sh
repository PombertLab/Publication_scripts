#!/usr/bin/bash

## fastp version 0.20.1

OUTDIR=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP
mkdir -p $OUTDIR

## 151 bp long
for x in {50451,50506,50507,50602,50603,50604,50651}; do 
	printf "Fastp working on $x...\n";
    fastp \
      -w 10 \
      -i $x.raw.R1.fastq.gz \
      -I $x.raw.R2.fastq.gz \
      -o $OUTDIR/$x.q30.l125.R1.fastq.gz \
      -O $OUTDIR/$x.q30.l125.R2.fastq.gz \
      -M 30 \
      -r \
      -l 125 \
	  --json $OUTDIR/$x.fastp.json \
   	  --html $OUTDIR/$x.fastp.html
done

## 125 bp long
printf "Fastp working on ECIII-L...\n"
fastp \
   -w 10 \
   -i ECIII_L.raw.R1.fastq.gz \
   -I ECIII_L.raw.R2.fastq.gz \
   -o $OUTDIR/ECIII_L.q30.l75.R1.fastq.gz \
   -O $OUTDIR/ECIII_L.q30.l75.R2.fastq.gz \
   -M 30 \
   -r \
   -l 75 \
   --json $OUTDIR/ECIII_L.fastp.json \
   --html $OUTDIR/ECIII_L.fastp.html

# 100 bp long
for x in {ECI,ECII,ECIII}; do 
	printf "Fastp working on $x...\n";
    fastp \
      -w 10 \
      -i $x.raw.R1.fastq.gz \
      -I $x.raw.R2.fastq.gz \
      -o $OUTDIR/$x.q30.l50.R1.fastq.gz \
      -O $OUTDIR/$x.q30.l50.R2.fastq.gz \
      -M 30 \
      -r \
      -l 50 \
	  --json $OUTDIR/$x.fastp.json \
   	  --html $OUTDIR/$x.fastp.html
done

## Swiss 54 bp long
printf "Fastp working on Swiss...\n"
fastp \
   -w 10 \
   -i Swiss.raw.R1.fastq.gz \
   -I Swiss.raw.R2.fastq.gz \
   -o $OUTDIR/Swiss.q30.l40.R1.fastq.gz \
   -O $OUTDIR/Swiss.q30.l40.R2.fastq.gz \
   -M 30 \
   -r \
   -l 40 \
   --json $OUTDIR/Swiss.fastp.json \
   --html $OUTDIR/Swiss.fastp.html

## ECII-CZ 50 bp long. single ends
printf "Fastp working on ECII-CZ...\n"
fastp \
   -w 10 \
   -i ECII_CZ.raw.SE.fastq.gz \
   -o $OUTDIR/ECII_CZ.q30.l40.SE.fastq.gz \
   -M 30 \
   -r \
   -l 40 \
   --json $OUTDIR/Swiss.fastp.json \
   --html $OUTDIR/Swiss.fastp.html

## Ordospora colligata; OC4 not available, 125 bp + 1
for x in {GBEP,NOV7,FISK}; do 
	printf "Fastp working on $x...\n";
    fastp \
      -w 10 \
      -i $x.raw.R1.fastq.gz \
      -I $x.raw.R2.fastq.gz \
      -o $OUTDIR/$x.q30.l125.R1.fastq.gz \
      -O $OUTDIR/$x.q30.l125.R2.fastq.gz \
      -M 30 \
      -r \
      -l 125 \
	  --json $OUTDIR/$x.fastp.json \
   	  --html $OUTDIR/$x.fastp.html
done
