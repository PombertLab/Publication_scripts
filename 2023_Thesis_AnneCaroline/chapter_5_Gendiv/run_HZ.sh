#!/usr/bin/env bash 

GENDIV=/home/jpombert/Analyses/Encephalitozoon/GenDiv/

## Setting --min-var-freq 0.2 to detect possible heterozigous diploid (var-freq ~ 0.5), 
## triploid (~ 0.33/0.66) or tetraploid (~ 0.25/0.5/0.75) sites

# Paired ends
for x in {50451,50604,Swiss,50506,50507,50603,50651,50602,ECIII_L,ECI,ECII,ECIII}; do 
	get_SNPs.pl \
	  -fa $GENDIV/FASTA/$x.fasta \
	  -pe1 $GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/$x.*.R1.fastq.gz \
	  -pe2 $GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/$x.*.R2.fastq.gz \
	  -mapper minimap2 \
	  -preset sr \
	  -caller varscan2 \
	  -var $GENDIV/VarScan.v2.4.6.jar \
	  -type both \
	  --min-var-freq 0.2 \
	  --min-avg-qual 30 \
	  -threads 12 \
	  -mem 2 \
	  -o $GENDIV/RM/HZ
done

# Paired ends Ordospora (OC4 not in SRA)
for x in {GBEP,NOV7,FISK}; do 
	get_SNPs.pl \
	  -fa $GENDIV/FASTA/$x.fasta \
	  -pe1 $GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/$x.*.R1.fastq.gz \
	  -pe2 $GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/$x.*.R2.fastq.gz \
	  -mapper minimap2 \
	  -preset sr \
	  -caller varscan2 \
	  -var $GENDIV/VarScan.v2.4.6.jar \
	  -type both \
	  --min-var-freq 0.2 \
	  --min-avg-qual 30 \
	  -threads 12 \
	  -mem 2 \
	  -o $GENDIV/RM/HZ
done

# Single ends (ECII-CZ)
get_SNPs.pl \
  -fa $GENDIV/FASTA/ECII_CZ.fasta \
  -fq $GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP/ECII_CZ.*.SE.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-var-freq 0.2 \
  --min-avg-qual 30 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/HZ

######## copy + rename VCFs for SnpEff
DIR=$GENDIV/snpEFF/VCFs_HZ/
SNPEFF=$GENDIV/snpEFF

$GENDIV/Scripts/rename_vcf.pl \
  -v $GENDIV/RM/HZ/minimap2.varscan2.VCFs/*.vcf \
  -o $DIR/

cd $DIR;

#####  Run snpEff with csv stats
#cd No 50504, GB_M1 or OC4 SRA data, skipping them...
for x in {50451,50604,Swiss,50506,50507,50603,50651,50602,ECII_CZ,ECIII_L,ECI,ECII,ECIII,NOV7,GBEP,FISK}; do
	REF=$x;
	VCF="$x.vs.$x.vcf";
	SUM="$x.vs.$x.summary.csv";
	OUT="$x.vs.$x.annotations.vcf";
	java -Xmx8g -jar $SNPEFF/snpEff.jar \
	  -c $SNPEFF/snpEff.config \
	  -v \
	  $REF \
	  $VCF \
	  -csvStats $SUM  \
	  > $OUT
done;

cd $GENDIV

## making the table
$GENDIV/Scripts/make_table_HZ.pl \
  -l $GENDIV/Scripts/Table_list_HZ.txt \
  -acc $GENDIV/Scripts/accessions.tsv \
  -stats $GENDIV/RM/HZ/minimap2.varscan2.stats/*.stats \
  -eff $GENDIV/snpEFF/VCFs_HZ/*.csv \
  -out $GENDIV/table_HZ.tsv
