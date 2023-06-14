#!/usr/bin/env bash

## Quick check for congruency with varscan2 results; should
## get similar number of SNPs... 

picard=/opt/gatk-4.2.6.1/picard.jar

for x in {50451,50504,50506,50507,50602,50604,50651,ECI,ECII,ECII_CZ,ECIII,ECIII_L,GB_M1,FISK,OC4,NOV7,GBEP,Swiss}; do
	java -jar $picard \
	  CreateSequenceDictionary \
	  -R $GENDIV/FASTA/$x.fasta \
	  -O $GENDIV/FASTA/$x.dict
done

get_SNPs.pl \
  -fa $GENDIV/FASTA/{50451,50504,50604,Swiss}.fasta \
  -pe1 $GENDIV/FASTQ/FASTQ_250/{50451,50504,50604,Swiss}.250.R1.fastq.gz \
  -pe2 $GENDIV/FASTQ/FASTQ_250/{50451,50504,50604,Swiss}.250.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -rmo \
  -bam \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/GATK_250

### Note haplotype caller realigns reads. No option to turn it off
# ploidy 2
for x in {50451,50504,50604,Swiss}; do
	for y in {50451,50504,50604,Swiss}; do
		gatk \
		  --java-options "-Xmx8G" \
		  HaplotypeCaller \
		  -ploidy 2 \
		  -R $GENDIV/FASTA/$x.fasta \
		  -I $GENDIV/RM/GATK_250/minimap2.BAM/$y.250.R1.fastq.gz.$x.fasta.minimap2.bam \
		  -O $GENDIV/RM/GATK_250/$y.vs.$x.gatk.vcf
	done
done

# ploidy 1
for x in {50451,50504,50604,Swiss}; do
	for y in {50451,50504,50604,Swiss}; do
		gatk \
		  --java-options "-Xmx8G" \
		  HaplotypeCaller \
		  -ploidy 1 \
		  -R $GENDIV/FASTA/$x.fasta \
		  -I $GENDIV/RM/GATK_250/minimap2.BAM/$y.250.R1.fastq.gz.$x.fasta.minimap2.bam \
		  -O $GENDIV/RM/GATK_250/$y.vs.$x.gatk.p1.vcf
	done
done

# ploidy 1 + no softclip
for x in {50451,50504,50604,Swiss}; do
	for y in {50451,50504,50604,Swiss}; do
		gatk \
		  --java-options "-Xmx8G" \
		  HaplotypeCaller \
		  -ploidy 1 \
		  --dont-use-soft-clipped-bases true \
		  -R $GENDIV/FASTA/$x.fasta \
		  -I $GENDIV/RM/GATK_250/minimap2.BAM/$y.250.R1.fastq.gz.$x.fasta.minimap2.bam \
		  -O $GENDIV/RM/GATK_250/$y.vs.$x.gatk.p1.nosoft.vcf
	done
done

