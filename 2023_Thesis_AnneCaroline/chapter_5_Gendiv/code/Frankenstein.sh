#!/usr/bin/env bash

##########################################################################################

# This test is to assess the effect of contamination from 2 distinct ppopulations
# on the outcome of the SNPs, heterozigosity and GCeq analyses

# Using 10/90, 20/80 and 30/70 split from similar (~1 SNP/kb) isolates
# 10M PE reads total

# Will help assess what is happening with the E. hellem Swiss isolate

##########################################################################################
### Creating chimeric datasets

FASTP=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP
FRK=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/FRANKENSTEIN

# k = number of reads (in millions)
for k in {1..4}; do 
	$GENDIV/Scripts/trim_fq.pl \
	  -f $FASTP/50602.q30.l125.R* \
	  -l 54 \
	  -n $k \
	  -o $FRK/DATA
done

for k in {6..9}; do 
	$GENDIV/Scripts/trim_fq.pl \
	  -f $FASTP/ECIII_L.q30.l75.R* \
	  -l 54 \
	  -n $k \
	  -o $FRK/DATA
done

cd $FRK/DATA

## FRK1 => 10% 50602 + 90% ECIII_L
cat $FRK/DATA/50602.q30.l125.R1.trimmed.1M.fastq  $FRK/DATA/ECIII_L.q30.l75.R1.trimmed.9M.fastq > $FRK/DATA/frk1.R1.fastq
cat 50602.q30.l125.R2.trimmed.1M.fastq  $FRK/DATA/ECIII_L.q30.l75.R2.trimmed.9M.fastq > $FRK/DATA/frk1.R2.fastq

## FRK2 => 20% 50602 + 80% ECIII_L
cat $FRK/DATA/50602.q30.l125.R1.trimmed.2M.fastq  $FRK/DATA/ECIII_L.q30.l75.R1.trimmed.8M.fastq > $FRK/DATA/frk2.R1.fastq
cat $FRK/DATA/50602.q30.l125.R2.trimmed.2M.fastq  $FRK/DATA/ECIII_L.q30.l75.R2.trimmed.8M.fastq > $FRK/DATA/frk2.R2.fastq

## FRK3 => 30% 50602 + 70% ECIII_L
cat $FRK/DATA/50602.q30.l125.R1.trimmed.3M.fastq  $FRK/DATA/ECIII_L.q30.l75.R1.trimmed.7M.fastq > $FRK/DATA/ frk3.R1.fastq
cat $FRK/DATA/50602.q30.l125.R2.trimmed.3M.fastq  $FRK/DATA/ECIII_L.q30.l75.R2.trimmed.7M.fastq > $FRK/DATA/frk3.R2.fastq

## FRK4 => 40% 50602 + 60% ECIII_L
cat $FRK/DATA/50602.q30.l125.R1.trimmed.4M.fastq  $FRK/DATA/ECIII_L.q30.l75.R1.trimmed.6M.fastq > $FRK/DATA/frk4.R1.fastq
cat $FRK/DATA/50602.q30.l125.R2.trimmed.4M.fastq  $FRK/DATA/ECIII_L.q30.l75.R2.trimmed.6M.fastq > $FRK/DATA/frk4.R2.fastq

## Clean up tmp files
rm $FRK/DATA/50602*
rm $FRK/DATA/ECIII*

##########################################################################################
### Assembling the chimeric datasets with SPAdes

mkdir -p $GENDIV/SPAdes

spades.py \
  -1 $FRK/DATA/frk1.R1.fastq \
  -2 $FRK/DATA/frk1.R2.fastq \
  --threads 12 \
  --memory 12 \
  --careful \
  -o $GENDIV/SPAdes/Frk1

spades.py \
  -1 $FRK/DATA/frk2.R1.fastq \
  -2 $FRK/DATA/frk2.R2.fastq \
  --threads 12 \
  --memory 12 \
  --careful \
  -o $GENDIV/SPAdes/Frk2

spades.py \
  -1 $FRK/DATA/frk3.R1.fastq \
  -2 $FRK/DATA/frk3.R2.fastq \
  --threads 12 \
  --memory 12 \
  --careful \
  -o $GENDIV/SPAdes/Frk3

spades.py \
  -1 $FRK/DATA/frk4.R1.fastq \
  -2 $FRK/DATA/frk4.R2.fastq \
  --threads 12 \
  --memory 12 \
  --careful \
  -o $GENDIV/SPAdes/Frk4

mkdir $GENDIV/FRANKEN

## Filtering out small contigs...
for x in {1..4}; do
	$GENDIV/Scripts/parse_FASTA.pl \
	 -f $GENDIV/SPAdes/Frk$x/contigs.fasta -m 5000 -o $GENDIV/FRANKEN/;
	mv $GENDIV/FRANKEN/contigs.5000.fasta $GENDIV/FRANKEN/frk$x.fasta
done

##########################################################################################
### Checking for contamination in assemblies; many low depth contigs according to spades

$GENDIV/Scripts/orient_fasta_to_reference.py \
  -f $GENDIV/FRANKEN/*.fasta \
  -r $GENDIV/FASTA/50602.fasta

# first hit in frk1.unmatched.fasta = mycoplasma! Very AT rich !

# Checking with Swiss E. hellem assembly too just in case:
$GENDIV/Scripts/orient_fasta_to_reference.py \
  -f $GENDIV/FASTA/Swiss.fasta \
  -r $GENDIV/FASTA/50451.fasta

# Not a problem (no contamination) so that is not the issue...

##########################################################################################
### Keeping only contigs matching to cuniculi (discarding contaminants)

for x in {1..4}; do
  cp $GENDIV/FRANKEN/oriented_fastas/frk$x/frk$x.oriented.fasta $GENDIV/FRANKEN/
done

##########################################################################################
### Read mapping with contaminated FASTQ datasets

mkdir -p $GENDIV/FRANKEN/RM
get_SNPs.pl \
  -fa $GENDIV/FASTA/{50602,GB_M1,ECII_CZ,ECIII_L,ECI,ECII,ECIII}.fasta \
  -pe1 $FRK/DATA/*.R1.fastq \
  -pe2 $FRK/DATA/*.R2.fastq \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.7 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/FRANKEN/RM/

# Read mapping HZ
mkdir -p $GENDIV/FRANKEN/HZ
get_SNPs.pl \
  -fa $GENDIV/FRANKEN/*.oriented.fasta \
  -pe1 $FRK/DATA/*.R1.fastq \
  -pe2 $FRK/DATA/*.R2.fastq \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  -bam \
  --min-var-freq 0.2 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/FRANKEN/HZ/

##########################################################################################
###  GCeq

GENDIV=/home/jpombert/Analyses/Encephalitozoon/GenDiv
VCF=$GENDIV/FRANKEN/RM/minimap2.varscan2.VCFs

$GENDIV/Scripts/gceq.pl \
  -f $GENDIV/FASTA/{50602,GB_M1,ECI,ECII,ECII_CZ,ECIII,ECIII_L}.fasta \
  -f $GENDIV/FRANKEN/*.oriented.fasta \
  -v \
  $VCF/*.vcf \
  -o $GENDIV/GCeq/Frk.gceq \
  -verbose


##########################################################################################
#### Ploidy NGS
mkdir -p $GENDIV/ploidyNGS

# Seems to stuggle with relative paths
cd /opt/ploidyNGS/

# Running them indepedently to use parallel threads
./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/Frk1_illumina \
  -b $GENDIV/FRANKEN/HZ/minimap2.BAM/frk1.R1.fastq.frk1.oriented.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/Frk2_illumina \
  -b $GENDIV/FRANKEN/HZ/minimap2.BAM/frk2.R1.fastq.frk2.oriented.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/Frk3_illumina \
  -b $GENDIV/FRANKEN/HZ/minimap2.BAM/frk3.R1.fastq.frk3.oriented.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/Frk4_illumina \
  -b $GENDIV/FRANKEN/HZ/minimap2.BAM/frk4.R1.fastq.frk4.oriented.fasta.minimap2.bam


##########################################################################################
#### R ploidy plots

VCF=/home/jpombert/Analyses/Encephalitozoon/GenDiv/FRANKEN/HZ/minimap2.varscan2.VCFs

$GENDIV/Scripts/sort_SNPs.pl \
  -vcf $VCF/*.vcf \
   -o $GENDIV/FRANKEN/Rplots \
   -noindel

$GENDIV/Scripts/plot_r.pl \
  -t $GENDIV/FRANKEN/Rplots/*.tsv \
  -o $GENDIV/FRANKEN/Rplots_PDFs \
  -ymax 0.12
