#!/usr/bin/env bash

############# Download nanopore data
mkdir -p $GENDIV/Nanopore
cd $GENDIV/Nanopore

## NCBI SRA accessions
# SRR23560420	Encephalitozoon hellem ATCC 50451
# SRR17853474	Encephalitozoon hellem ATCC 50604
# SRR17865590	Encephalitozoon intestinalis ATCC 50506
# SRR17858634	Encephalitozoon cuniculi ATCC 50602

for x in {SRR23560420,SRR17853474,SRR17865590,SRR17858634}; do
	fasterq-dump $x
done

## Gzip compression
pigz *.fastq

## Checking quality scores with FastQC
fastqc -t 8 $GENDIV/Nanopore/*.fastq.gz

############# Mapping with minimap2
## Calling HZ variants to see what happens; good or bad?
## Average qual estimated from FastQC: ~ 18 or so
## mpileup was choking due to the BAQ alignment of nanopore reads
## had to implement a noBAQ cmd line switch in get_SNPs.pl to prevent hiccups

# 50451; SRR23560420; Encephalitozoon hellem ATCC 50451
get_SNPs.pl \
  -fa $GENDIV/FASTA/50451.fasta \
  -fq SRR23560420.fastq.gz \
  -mapper minimap2 \
  -preset map-ont \
  -bam \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-avg-qual 12 \
  --min-var-freq 0.2 \
  --noBAQ \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/Nanopore

# 50604; SRR17853474; Encephalitozoon hellem ATCC 50604
get_SNPs.pl \
  -fa $GENDIV/FASTA/50604.fasta \
  -fq SRR17853474.fastq.gz \
  -mapper minimap2 \
  -preset map-ont \
  -bam \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-avg-qual 12 \
  --min-var-freq 0.2 \
  --noBAQ \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/Nanopore

# 50506; SRR17865590; Encephalitozoon intestinalis ATCC 50506
get_SNPs.pl \
  -fa $GENDIV/FASTA/50506.fasta \
  -fq SRR17865590.fastq.gz \
  -mapper minimap2 \
  -preset map-ont \
  -bam \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-avg-qual 12 \
  --min-var-freq 0.2 \
  --noBAQ \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/Nanopore

# 50602; SRR17858634; Encephalitozoon cuniculi ATCC 50602
get_SNPs.pl \
  -fa $GENDIV/FASTA/50602.fasta \
  -fq SRR17858634.fastq.gz \
  -mapper minimap2 \
  -preset map-ont \
  -bam \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-avg-qual 12 \
  --min-var-freq 0.2 \
  --noBAQ \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/Nanopore

## Keeping only the primary alignments from the bam files
## Doesnt appear to affect the plots from
## ploidyNGS, same plots with or without them
## Doing this step doesnt improve things

mkdir -p $GENDIV/Nanopore/Primary_align/

# 50451
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  $GENDIV/Nanopore/minimap2.BAM/SRR23560420.fastq.gz.50451.fasta.minimap2.bam \
  > $GENDIV/Nanopore/Primary_align/50451.primary.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/Nanopore/Primary_align/50451.primary.sorted.bam \
  $GENDIV/Nanopore/Primary_align/50451.primary.bam

# 50604
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  $GENDIV/Nanopore/minimap2.BAM/SRR17853474.fastq.gz.50604.fasta.minimap2.bam \
  > $GENDIV/Nanopore/Primary_align/50604.primary.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/Nanopore/Primary_align/50604.primary.sorted.bam \
  $GENDIV/Nanopore/Primary_align/50604.primary.bam

# 50506
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  $GENDIV/Nanopore/minimap2.BAM/SRR17865590.fastq.gz.50506.fasta.minimap2.bam \
  > $GENDIV/Nanopore/Primary_align/50506.primary.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/Nanopore/Primary_align/50506.primary.sorted.bam \
  $GENDIV/Nanopore/Primary_align/50506.primary.bam

# 50602
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  -b $GENDIV/Nanopore/minimap2.BAM/SRR17858634.fastq.gz.50602.fasta.minimap2.bam \
   > $GENDIV/Nanopore/Primary_align/50602.primary.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/Nanopore/Primary_align/50602.primary.sorted.bam \
  $GENDIV/Nanopore/Primary_align/50602.primary.bam

#### Ploidy NGS
mkdir -p $GENDIV/ploidyNGS

# Seems to stuggle with relative paths
cd /opt/ploidyNGS/

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50451_nanopore \
  -b $GENDIV/Nanopore/Primary_align/50451.primary.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50604_nanopore \
  -b $GENDIV/Nanopore/Primary_align/50604.primary.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50506_nanopore \
  -b $GENDIV/Nanopore/Primary_align/50506.primary.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50602_nanopore \
  -b $GENDIV/Nanopore/Primary_align/50602.primary.sorted.bam



### Same test with illumina (removing secondary alignments)
### The raw BAM files is messy, would doing this get rid of 
### background noise? Answer: no
#
### Removing secondary alignments has no effect on
### ploidyNGS and varscan outputs, those are already ignored
### by the tools

HZ=$GENDIV/RM/HZ
mkdir -p $HZ/Primary_align
mkdir -p $GENDIV/RM/HZ/CC

## 50451
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  $HZ/minimap2.BAM/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.bam \
  > $HZ/Primary_align/50451.primary.bam

samtools sort \
  -@ 8 \
  -o $HZ/Primary_align/50451.primary.sorted.bam \
  $HZ/Primary_align/50451.primary.bam

## Swiss
samtools view \
  -bS \
  -F 256 \
  -@ 8 \
  $HZ/minimap2.BAM/Swiss.q30.l40.R1.fastq.gz.Swiss.fasta.minimap2.bam \
  > $HZ/Primary_align/Swiss.primary.bam

samtools sort \
  -@ 8 \
  -o $HZ/Primary_align/Swiss.primary.sorted.bam \
  $HZ/Primary_align/Swiss.primary.bam

###########################   BAM CC subsets   ##############################
# 50451
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.bam \
	"CP119062.1:25000-235000" \
	"CP119063.1:25000-200000" \
	"CP119064.1:20000-200000" \
	"CP119065.1:30000-215000" \
	"CP119066.1:30000-220000" \
	"CP119067.1:30000-225000" \
	"CP119068.1:30000-230000" \
	"CP119069.1:20000-210000" \
	"CP119070.1:20000-185000" \
	"CP119071.1:30000-250000" \
	"CP119072.1:20000-245000" \
	> $GENDIV/RM/HZ/CC/50451.cc.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/RM/HZ/CC/50451.cc.sorted.bam \
  $GENDIV/RM/HZ/CC/50451.cc.bam

# 50604
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.bam \
	"CP075147.1:35000-240000" \
	"CP075148.1:30000-190000" \
	"CP075149.1:25000-195000" \
	"CP075150.1:35000-220000" \
	"CP075151.1:30000-220000" \
	"CP075152.1:30000-225000" \
	"CP075153.1:30000-230000" \
	"CP075154.1:30000-220000" \
	"CP075155.1:30000-190000" \
	"CP075156.1:27500-250000" \
	"CP075157.1:25000-245000" \
	> $GENDIV/RM/HZ/CC/50604.cc.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/RM/HZ/CC/50604.cc.sorted.bam \
  $GENDIV/RM/HZ/CC/50604.cc.bam

# 50506
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.bam \
	"CP075158.1:25000-160000" \
	"CP075159.1:30000-185000" \
	"CP075160.1:20000-190000" \
	"CP075161.1:30000-210000" \
	"CP075162.1:20000-210000" \
	"CP075163.1:20000-210000" \
	"CP075164.1:15000-220000" \
	"CP075165.1:30000-220000" \
	"CP075166.1:20000-250000" \
	"CP075167.1:30000-255000" \
	"CP075168.1:30000-250000" \
	> $GENDIV/RM/HZ/CC/50506.cc.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/RM/HZ/CC/50506.cc.sorted.bam \
  $GENDIV/RM/HZ/CC/50506.cc.bam

# 50602
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.bam \
	"CP091441.1:40000-180000" \
	"CP091440.1:40000-200000" \
	"CP091439.1:30000-210000" \
	"CP091436.1:30000-220000" \
	"CP091438.1:25000-210000" \
	"CP091437.1:20000-215000" \
	"CP091435.1:25000-225000" \
	"CP091434.1:40000-240000" \
	"CP091433.1:25000-265000" \
	"CP091432.1:30000-270000" \
	"CP091431.1:50000-270000" \
	> $GENDIV/RM/HZ/CC/50602.cc.bam

samtools sort \
  -@ 8 \
  -o $GENDIV/RM/HZ/CC/50602.cc.sorted.bam \
  $GENDIV/RM/HZ/CC/50602.cc.bam


##############################   ploidyNGS   ################################
cd /opt/ploidyNGS/

## Raw
./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50506_illumina \
  -b $HZ/minimap2.BAM/50506.q30.l125.R1.fastq.gz.50506.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50602_illumina \
  -b $HZ/minimap2.BAM/50602.q30.l125.R1.fastq.gz.50602.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50604_illumina \
  -b $HZ/minimap2.BAM/50604.q30.l125.R1.fastq.gz.50604.fasta.minimap2.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/Swiss_illumina \
  -b $HZ/Primary_align/Swiss.primary.sorted.bam

## Test with primary alignments only, no change observed
./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50451_ilumina \
  -b $HZ/Primary_align/50451.primary.sorted.bam

## CC subsets
./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50451_illumina_cc \
  -b $GENDIV/RM/HZ/CC/50451.cc.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50506_illumina_cc \
  -b $GENDIV/RM/HZ/CC/50506.cc.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50602_illumina_cc \
  -b $GENDIV/RM/HZ/CC/50602.cc.sorted.bam

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50604_illumina_cc \
  -b $GENDIV/RM/HZ/CC/50604.cc.sorted.bam


### Test on data mapped on FASTA deaturing only the CC subset
### Not great, just maps repeats at other locations??
FASTP=$GENDIV/FASTQ/FASTQ_ILLUMINA/FASTP

get_SNPs.pl \
  -fa $GENDIV/FASTA/CC/50451.cc.fasta \
  -pe1 $FASTP/50451.q30.l125.R1.fastq.gz \
  -pe2 $FASTP/50451.q30.l125.R2.fastq.gz \
  -mapper minimap2 \
  -preset sr \
  -bam \
  -caller varscan2 \
  -var $GENDIV/VarScan.v2.4.6.jar \
  -type both \
  --min-avg-qual 20 \
  --min-var-freq 0.2 \
  -threads 12 \
  -mem 2 \
  -o $GENDIV/RM/CC

./ploidyNGS.py \
  -o $GENDIV/ploidyNGS/50451_ilumina_cc \
  -b $GENDIV/RM/CC/minimap2.BAM/50451.q30.l125.R1.fastq.gz.50451.cc.fasta.minimap2.bam







### 
 2010  mv SRR23560420.fastq.gz.50451.fasta.minimap2.both.vcf test2.vcf
 2011  mv test2.vcf 50451.vcf
 $GENDIV/Scripts/parse_VCF_by_CC.pl -l $GENDIV/Scripts/cc_list.txt -v 50451.vcf -o CC


## check rRNA and PTP1 50451
$GENDIV/Scripts/extract_subfasta.pl \
  -f $GENDIV/FASTA/50451.fasta \
  -l \
  CP119062.1:3075-7173 \
  CP119067.1:51565-52986 \
  -o $GENDIV/TEST/BAM_subsets

# 50451 rRNA
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.bam \
	"CP119062.1:3075-7173" \
	> $GENDIV/TEST/50451.rRNA.bam

samtools index $GENDIV/TEST/BAM_subsets/50451.rRNA.bam

# 50451 PTP1
samtools view \
	-bS \
	-@ 8 \
	$GENDIV/RM/HZ/minimap2.BAM/50451.q30.l125.R1.fastq.gz.50451.fasta.minimap2.bam \
	"CP119067.1:51565-52986" \
	> $GENDIV/TEST/BAM_subsets/50451.PTP1.bam

samtools index $GENDIV/TEST/50451.PTP1.bam