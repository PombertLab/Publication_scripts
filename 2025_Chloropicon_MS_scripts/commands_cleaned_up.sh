#!/usr/bin/env bash

####################################################################################
### Creating a project folder
####################################################################################

CP=/mnt/d/LinuxVM/Chloropicon/

mkdir -p $CP
mkdir -p $CP/FASTA
mkdir -p $CP/READS
mkdir -p $CP/READS/RAW
mkdir -p $CP/READS/RAW/CCMP1205   ## Data from NCBI SRA (Lemieux et al. DOI: 10.1038/s41467-019-12014-x)
mkdir -p $CP/READS/RAW/RCC138     ## Data deposited in NCBI SRA (this study)
mkdir -p $CP/READS/RAW/RCC2335    ## Data deposited in NCBI SRA (this study)
mkdir -p $CP/READS/RAW/CCMP1998   ## Data deposited in NCBI SRA (this study)
mkdir -p $CP/READS/RAW/RCC1871    ## Data deposited in NCBI SRA (this study)
mkdir -p $CP/STRUCTURES           ## For structural homology
mkdir -p $CP/Scripts


####################################################################################
### Sequencing data
####################################################################################

######## NCBI SRA toolkit:
# https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

############### CCMP1205 ###############

### Downloading Chloropicon primus CCMP1205 NCBI SRA accessions:
# SRR8185497  Illumina PE
# SRR8185495  PacBio # 1
# SRR8185494  PacBio # 2
# SRR8185493  PacBio # 3
# SRR8185492  PacBio # 4

for x in {SRR8185497,SRR8185495,SRR8185494,SRR8185493,SRR8185492}; do
  echo "Downloading dataset $x from the NCBI SRA archive..."
  fasterq-dump $x \
  --outdir $CP/READS/RAW/CCMP1205 \
  --progress
done

# Removing tmp file left, if any
rm -R $CP/READS/RAW/CCMP1205/fasterq.tmp*

# Illumina
mv $CP/READS/RAW/CCMP1205/SRR8185497_1.fastq $CP/READS/RAW/CCMP1205/CCMP1205_PE_R1.fastq
mv $CP/READS/RAW/CCMP1205/SRR8185497_2.fastq $CP/READS/RAW/CCMP1205/CCMP1205_PE_R2.fastq

pigz $CP/READS/RAW/CCMP1205/CCMP1205_R1.fastq
pigz $CP/READS/RAW/CCMP1205/CCMP1205_R2.fastq

# PacBio
cat \
  $CP/READS/RAW/CCMP1205/SRR8185495.fastq \
  $CP/READS/RAW/CCMP1205/SRR8185494.fastq \
  $CP/READS/RAW/CCMP1205/SRR8185493.fastq \
  $CP/READS/RAW/CCMP1205/SRR8185492.fastq \
  > $CP/READS/RAW/CCMP1205/CCMP1205.pacbio.fastq

pigz $CP/READS/RAW/CCMP1205/CCMP1205.pacbio.fastq

# Cleanup
rm $CP/READS/RAW/CCMP1205/SRR*.fastq


############### RCC1871 ###############
for x in {SRR27893523,SRR27893524,SRR27893525,SRR27893526}; do
  echo "Downloading dataset $x from the NCBI SRA archive..."
  fasterq-dump $x \
    --outdir $CP/READS/RAW/RCC1871 \
    --progress
done

# PacBio
cat \
  $CP/READS/RAW/RCC1871/SRR27893523.fastq \
  $CP/READS/RAW/RCC1871/SRR27893524.fastq \
  > $CP/READS/RAW/RCC1871/RCC1871.pacbio.fastq

pigz $CP/READS/RAW/RCC1871/RCC1871.pacbio.fastq

# Illumina
for x in {1,2}; do
  cat $CP/READS/RAW/RCC1871/SRR27893525_$x.fastq \
    $CP/READS/RAW/RCC1871/SRR27893526_$x.fastq \
    > $CP/READS/RAW/RCC1871/RCC1871_R$x.fastq
  pigz $CP/READS/RAW/RCC1871/RCC1871_R$x.fastq
done

# Cleanup
rm $CP/READS/RAW/RCC1871/SRR*.fastq

############### CCMP1998 ###############
for x in {SRR27907038,SRR27907039,SRR27907040}; do
  echo "Downloading dataset $x from the NCBI SRA archive..."
  fasterq-dump $x \
    --outdir $CP/READS/RAW/CCMP1998 \
    --progress
done

# Pacbio
cat \
  $CP/READS/RAW/CCMP1998/SRR27907038.fastq \
  $CP/READS/RAW/CCMP1998/SRR27907039.fastq \
  > $CP/READS/RAW/CCMP1998/CCMP1998.pacbio.fastq

pigz $CP/READS/RAW/CCMP1998/CCMP1998.pacbio.fastq

# Illumina
for x in {1,2}; do
  mv $CP/READS/RAW/CCMP1998/SRR27907040_$x.fastq $CP/READS/RAW/CCMP1998/CCMP1998_R$x.fastq
  pigz $CP/READS/RAW/CCMP1998/CCMP1998_R$x.fastq
done

# Cleanup
rm $CP/READS/RAW/CCMP1998/SRR27907038.fastq

############### RCC2335 ###############
for x in {SRR27895191,SRR27895192,SRR27895193,SRR27895194}; do
  echo "Downloading dataset $x from the NCBI SRA archive..."
  fasterq-dump $x \
    --outdir $CP/READS/RAW/RCC2335 \
    --progress
done

# PacBio/Nanopore
cat \
  $CP/READS/RAW/RCC2335/SRR27895191.fastq \
  $CP/READS/RAW/RCC2335/SRR27895192.fastq \
  > $CP/READS/RAW/RCC2335/RCC2335.pacbio.fastq

mv $CP/READS/RAW/RCC2335/SRR27895193.fastq $CP/READS/RAW/RCC2335/RCC2335.nanopore.fastq

pigz $CP/READS/RAW/RCC2335/RCC2335.pacbio.fastq
pigz $CP/READS/RAW/RCC2335/RCC2335.nanopore.fastq

# Illumina
for x in {1,2}; do
  mv $CP/READS/RAW/RCC2335/SRR27895194_$x.fastq $CP/READS/RAW/RCC2335/RCC2335_R$x.fastq
  pigz $CP/READS/RAW/RCC2335/RCC2335_R$x.fastq
done

# Cleanup
rm $CP/READS/RAW/RCC2335/SRR*.fastq

############### RCC138 ###############
for x in {SRR27893602,SRR27893603,SRR27893604}; do
  echo "Downloading dataset $x from the NCBI SRA archive..."
  fasterq-dump $x \
    --outdir $CP/READS/RAW/RCC138 \
    --progress
done

# Nanopore
mv $CP/READS/RAW/RCC138/SRR27893602.fastq $CP/READS/RAW/RCC138/RCC138.nanopore.fastq
pigz $CP/READS/RAW/RCC138/RCC138.nanopore.fastq

# Illumina
for x in {1,2}; do
  cat \
    $CP/READS/RAW/RCC138/SRR27893603_$x.fastq \
    $CP/READS/RAW/RCC138/SRR27893604_$x.fastq \
    > $CP/READS/RAW/RCC138/RCC138_R$x.fastq
  pigz $CP/READS/RAW/RCC138/RCC138_R$x.fastq
done

# Cleanup
rm $CP/READS/RAW/RCC138/SRR*.fastq


############### RNAseq datasets ###############
# CCMP1998
fasterq-dump SRR1300422 \
  --outdir $CP/READS/RAW/CCMP1998 \
  --progress

mv $CP/READS/RAW/CCMP1998/SRR1300422_1.fastq $CP/READS/RAW/CCMP1998/CCMP1998_R1.RNAseq.fastq
mv $CP/READS/RAW/CCMP1998/SRR1300422_2.fastq $CP/READS/RAW/CCMP1998/CCMP1998_R2.RNAseq.fastq

pigz $CP/READS/RAW/CCMP1998/CCMP1998_R1.RNAseq.fastq
pigz $CP/READS/RAW/CCMP1998/CCMP1998_R2.RNAseq.fastq

# RCC1871
fasterq-dump SRR27893522 \
  --outdir $CP/READS/RAW/RCC1871 \
  --progress

mv $CP/READS/RAW/RCC1871/SRR27893522_1.fastq $CP/READS/RAW/RCC1871/RCC1871_R1.RNAseq.fastq
mv $CP/READS/RAW/RCC1871/SRR27893522_2.fastq $CP/READS/RAW/RCC1871/RCC1871_R2.RNAseq.fastq

pigz $CP/READS/RAW/RCC1871/RCC1871_R1.RNAseq.fastq
pigz $CP/READS/RAW/RCC1871/RCC1871_R2.RNAseq.fastq

# RCC2335
fasterq-dump SRR27895190 \
  --outdir $CP/READS/RAW/RCC2335 \
  --progress

mv $CP/READS/RAW/RCC2335/SRR27895190_1.fastq $CP/READS/RAW/RCC2335/RCC2335_R1.RNAseq.fastq
mv $CP/READS/RAW/RCC2335/SRR27895190_2.fastq $CP/READS/RAW/RCC2335/RCC2335_R2.RNAseq.fastq

pigz $CP/READS/RAW/RCC2335/RCC2335_R1.RNAseq.fastq
pigz $CP/READS/RAW/RCC2335/RCC2335_R2.RNAseq.fastq

####################################################################################
### Gene colinearity/genome colinearity (synteny) analyses with SYNY
####################################################################################

mkdir -p $CP/SYNY

ls $CP/GBFF  # CCMP1205.gbff  CCMP1998.gbff  RCC138.gbff  RCC1871.gbff  RCC2335.gbff

### CCMP1205 vs. ALL
run_syny.pl \
  -a $CP/GBFF/*.gbff \
  -g 0 1 5 \
  -e 1e-10 \
  -o $CP/SYNY/ALL \
  -r CCMP1205 \
  --custom_preset chloropicon

### CCMP1205 vs. RCC138
run_syny.pl \
  -a $CP/GBFF/{CCMP1205,RCC138}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/CCMP1205_vs_RCC138 \
  -r CCMP1205 \
  --custom_preset chloropicon

### CCMP1205 vs. RCC2335
run_syny.pl \
  -a $CP/GBFF/{CCMP1205,RCC2335}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/CCMP_vs_RCC2335 \
  -r CCMP1205 \
  --custom_preset chloropicon

### CCMP1205 vs. RCC1871
run_syny.pl \
  -a $CP/GBFF/{CCMP1205,RCC1871}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/CCMP1205_vs_RCC1871 \
  -r CCMP1205 \
  --custom_preset chloropicon

### CCMP1205 vs. CCMP1998
run_syny.pl \
  -a $CP/GBFF/{CCMP1205,CCMP1998}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/CCMP1205_vs_CCMP1998 \
  -r CCMP1205 \
  --custom_preset chloropicon

### RCC2335 vs. CCMP1998
run_syny.pl \
  -a $CP/GBFF/{RCC2335,CCMP1998}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC2335_vs_CCMP1998 \
  -r RCC2335 \
  --custom_preset chloropicon

### RCC2335 vs. RCC1871
run_syny.pl \
  -a $CP/GBFF/{RCC2335,RCC1871}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC2335_vs_RCC1871 \
  -r RCC2335 \
  --custom_preset chloropicon

### RCC2335 vs. CCMP1205
run_syny.pl \
  -a $CP/GBFF/{RCC2335,CCMP1205}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC2335_vs_CCMP1205 \
  -r RCC2335 \
  --custom_preset chloropicon

### RCC2335 vs. RCC138
run_syny.pl \
  -a $CP/GBFF/{RCC2335,RCC138}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC2335_vs_RCC138 \
  -r RCC2335 \
  --custom_preset chloropicon

### RCC1871 vs. CCMP1998
run_syny.pl \
  -a $CP/GBFF/{RCC1871,CCMP1998}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC1871_vs_CCMP1998 \
  -r RCC1871 \
  --custom_preset chloropicon

### RCC1871 vs. RCC2335
run_syny.pl \
  -a $CP/GBFF/{RCC1871,RCC2335}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC1871_vs_RCC2335 \
  -r RCC1871 \
  --custom_preset chloropicon

### RCC1871 vs. CCMP1205
run_syny.pl \
  -a $CP/GBFF/{RCC1871,CCMP1205}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC1871_vs_CCMP1205 \
  -r RCC1871 \
  --custom_preset chloropicon

### RCC1871 vs. RCC138
run_syny.pl \
  -a $CP/GBFF/{RCC1871,RCC138}.gbff \
  -g 0 1 5 \
  -o $CP/SYNY/RCC1871_vs_RCC138 \
  -r RCC1871 \
  --custom_preset chloropicon


## Change the desired links file in concatenated.conf (if desired)
## => links file in concatenated.conf will default to the first one from -g option (above: .gap_0.links)
## Made some tweaks to the Circos configurations files to add ploidy (see below)
## then created final figures

####################################################################################
### fastANI 1.34
####################################################################################

run_fastANI.pl \
  -q $CP/FASTA/*.fasta \
  -r $CP/FASTA/*.fasta \
  -o $CP/FASTA//fastANI \
  -t 8

####################################################################################
### BUSCO 5.6.1
####################################################################################

mkdir -p $CP/BUSCO
mkdir -p $CP/BUSCO/Summaries
cd $CP/BUSCO

## Using protein sequences: SYNY/ALL/PROT_SEQ/

for file in {CCMP1205,CCMP1998,RCC138,RCC1871,RCC2335}; do
  busco \
    --mode proteins \
    --in $CP/SYNY/ALL/PROT_SEQ/$file.faa \
    --out $file \
    --out_path $CP/BUSCO/ \
    --lineage_dataset chlorophyta_odb10 \
    --cpu 8
done

## Copying summaries

for file in {CCMP1205,CCMP1998,RCC138,RCC1871,RCC2335}; do
  cp \
    $CP/BUSCO/$file/short_summary.*.txt \
    $CP/BUSCO/Summaries/
done

## Plotting summaries
BUSCOSC=/opt/busco/scripts/

python3 \
  $BUSCOSC/generate_plot.py \
  -wd $CP/BUSCO/Summaries/

####################################################################################
### Masking genomes (RepeatModeler 2.0.5 + RepeatMasker 4.1.5)
####################################################################################

mkdir -p $CP/REPEATMASKER

for species in {CCMP1998,CCMP1205,RCC138,RCC1871,RCC2335}; do

  mkdir -p $CP/REPEATMASKER/$species;
  cd $CP/REPEATMASKER/$species

  BuildDatabase \
    -name $species \
    $CP/FASTA/$species.fasta;

  echo "Running RepeatModeler on $species..."
  RepeatModeler \
    -database $species \
    -LTRStruct \
    -threads 8 \
    > rmblast.$species.out;

  ## RepeatMasker chokes if using hmmer as engine with fasta 
  ## file as lib input...
  echo "Running RepeatMasker on $species..."
  RepeatMasker \
    -e rmblast \
    -no_is \
    -lib ${species}-families.fa \
    $CP/FASTA/$species.fasta

done

mkdir -p $CP/FASTA/Masked

for species in {CCMP1998,CCMP1205,RCC138,RCC1871,RCC2335}; do
  mv \
    $CP/FASTA/$species.fasta.masked \
    $CP/FASTA/Masked/$species.masked.fasta
done

## Count masked regions (Nns) in sequences (Repeat % in Table 1)
$CP/Scripts/count_biases.pl \
  -f $CP/FASTA/Masked/*.fasta \
  -o $CP/FASTA/Masked/masked.per_contig.tsv \
  -s $CP/FASTA/Masked/masked.summary.tsv

## Metrics for Circos
for species in {CCMP1998,CCMP1205,RCC138,RCC1871,RCC2335}; do
  ~/Github/SYNY/Circos/nucleotide_biases.pl \
    -fasta $CP/FASTA/masked/$species.masked.fasta \
    -o $CP/REPEATMASKER/metrics \
    -winsize 10000 \
    -step 5000 \
    -ncheck \
    -custom
done

cat \
  $CP/REPEATMASKER/metrics/CCMP1205.masked/CCMP1205.masked.NN \
  $CP/REPEATMASKER/metrics/CCMP1998.masked/CCMP1998.masked.NN \
  $CP/REPEATMASKER/metrics/RCC138.masked/RCC138.masked.NN \
  $CP/REPEATMASKER/metrics/RCC1871.masked/RCC1871.masked.NN \
  $CP/REPEATMASKER/metrics/RCC2335.masked/RCC2335.masked.NN \
  > $CP/REPEATMASKER/metrics/concatenated.nn


####################################################################################
### Ploidy analysis (leveraging samtools, varscan2, R, and SSRG's get_SNPs.pl)
####################################################################################

### Assessing ploidy per chromosome by looking at the sequencing depths on each
### chromosome and by checking the allelic frequencies calculated with varscan2
### from the read mapping

mkdir -p $CP/PLOIDY
cd $CP/PLOIDY

### Copying genomes extracted with SYNY + getting rid of added RCC2335_ prefix
cp $CP/SYNY/CCMP1205_vs_RCC138/GENOME/*.fasta $CP/FASTA/
cp $CP/SYNY/CCMP1205_vs_RCC2335/GENOME/*.fasta $CP/FASTA/
sed 's/RCC2335_//' $CP/SYNY/CCMP_vs_RCC138/GENOME/RCC2335.fasta $CP/FASTA/RCC2335.fasta


####################################################################################
### Getting Illumina read quality metrics with FastQC v0.12.1
####################################################################################

mkdir -p $CP/READS/FastQC

## Genomic data
for species in {RCC138,RCC2335,CCMP1205,RCC1871,CCMP1998}; do 
  fastqc \
    --threads 8 \
    $CP/READS/RAW/$species/${species}_R1.fastq.gz \
    $CP/READS/RAW/$species/${species}_R2.fastq.gz \
    -o $CP/READS/FastQC
done

## RNAseq data
for species in {RCC1871,CCMP1998,RCC2335}; do
    fastqc \
    --threads 8 \
    $CP/READS/RAW/$species/${species}_R1.RNAseq.fastq.gz \
    $CP/READS/RAW/$species/${species}_R2.RNAseq.fastq.gz \
    -o $CP/READS/FastQC
done

####################################################################################
### Filtering the Illumina reads with FASTP 0.23.4
####################################################################################

mkdir -p $CP/READS/FASTP

## Genomic data
for species in {RCC138,RCC2335,CCMP1205,RCC1871,CCMP1998}; do 
  fastp \
    -w 10 \
    -i $CP/READS/RAW/$species/${species}_R1.fastq.gz \
    -I $CP/READS/RAW/$species/${species}_R2.fastq.gz \
    -o $CP/READS/FASTP/${species}_R1.q30.l125.fastq.gz \
    -O $CP/READS/FASTP/${species}_R2.q30.l125.fastq.gz \
    -M 30 \
    -r \
    -l 125 \
    --adapter_sequence=CTGTCTCTTATACACATCT \
    --adapter_sequence_r2=CTGTCTCTTATACACATCT \
    -j $CP/READS/FASTP/${species}_fastp.json \
    -h $CP/READS/FASTP/${species}_fastp.html 
done

## RNAseq data
for species in {RCC1871,CCMP1998,RCC2335}; do 
  fastp \
    -w 10 \
    -i $CP/READS/RAW/$species/${species}_R1.RNAseq.fastq.gz \
    -I $CP/READS/RAW/$species/${species}_R1.RNAseq.fastq.gz \
    -o $CP/READS/FASTP/${species}_R1.q30.l125.RNAseq.fastq.gz \
    -O $CP/READS/FASTP/${species}_R2.q30.l125.RNAseq.fastq.gz \
    -M 30 \
    -r \
    -l 125 \
    --adapter_sequence=CTGTCTCTTATACACATCT \
    --adapter_sequence_r2=CTGTCTCTTATACACATCT \
    -j $CP/READS/FASTP/${species}_RNAseq_fastp.json \
    -h $CP/READS/FASTP/${species}_RNAseq_fastp.html 
done

####################################################################################
### Creating data summary with MultiQC 1.19
####################################################################################

## pip install multiqc
## pip install QUAST
## pip install pragzip    ### for read_len_plot.py

mkdir -p $CP/MultiQC
mkdir -p $CP/MultiQC/BUSCO
mkdir -p $CP/MultiQC/FASTP
mkdir -p $CP/MultiQC/KLR
mkdir -p $CP/MultiQC/PLOTS

##### Quick QUAST report
quast.py \
  -o $CP/MultiQC/QUAST \
  $CP/FASTA/*.fasta

##### Copying BUSCO summaries
cp $CP/BUSCO/Summaries/* $CP/MultiQC/BUSCO/

##### Copying FASTQC reports
cp -R $CP/READS/FASTQC $CP/MultiQC/

##### Copying FASTP json/html reports
cp -R $CP/READS/FASTP/*.{json,html} $CP/MultiQC/FASTP

##### Long read metrics (_mqc.json file) for multiqc table with KLR
for species in {RCC1871,CCMP1205,CCMP1998}; do
  keep_longest_reads.pl \
    -input $CP/READS/RAW/$species/$species.pacbio.fastq.gz \
    -outdir $CP/MultiQC/KLR/$species \
    -metrics \
    -tsv \
    -json \
    -head 'Long read data'
done

keep_longest_reads.pl \
  -input $CP/READS/RAW/RCC138/RCC138.nanopore.fastq.gz \
  -outdir $CP/MultiQC/KLR/RCC138 \
  -metrics \
  -tsv \
  -json \
  -head 'Long read data'

keep_longest_reads.pl \
  -input $CP/READS/RAW/RCC2335/RCC2335.{pacbio,nanopore}.fastq.gz \
  -outdir $CP/MultiQC/KLR/RCC2335 \
  -metrics \
  -tsv \
  -json \
  -head 'Long read data'

##### PacBio/Nanopore read length distributions
for species in {CCMP1205,RCC1871,CCMP1998,RCC2335}; do
  read_len_plot.py \
    -f $CP/READS/RAW/${species}/${species}.pacbio.fastq.gz \
    --outdir  $CP/MultiQC/PLOTS \
    --output  ${species}_pacbio_read_distribution.svg \
              ${species}_pacbio_read_distribution.pdf \
              ${species}_pacbio_mqc.png \
    --metrics ${species}_pacbio_read_metrics.txt \
    -c blue \
    -x 50000 \
    --verbose
done

for species in {RCC138,RCC2335}; do
  read_len_plot.py \
    -f $CP/READS/RAW/${species}/${species}.nanopore.fastq.gz \
    --outdir  $CP/MultiQC/PLOTS \
    --output  ${species}_nanopore_read_distribution.svg \
              ${species}_nanopore_read_distribution.pdf \
              ${species}_nanopore_mqc.png \
    --metrics ${species}_nanopore_read_metrics.txt \
    -c blue \
    -x 50000 \
    --verbose
done

### Creating MultiQC HTML report summary
multiqc \
  --force \
  --outdir $CP/MultiQC/Summary \
  $CP/MultiQC

####################################################################################
### Performing the read mapping and variant calling step
####################################################################################

## Note: use a min var frequency of 0.1 => default [0.7] in
## get_SNPs.pl is for comparing bacterial genomes (haploid)

varscan=/opt/varscan2/VarScan.v2.4.6.jar

## Unmasked genomes
for species in {RCC138,RCC2335,RCC1871,CCMP1205,CCMP1998}; do 
  get_SNPs.pl \
    --fasta $CP/FASTA/$species.fasta \
    --pe1 $CP/READS/FASTP/${species}_R1.q30.l125.fastq.gz \
    --pe2 $CP/READS/FASTP/${species}_R2.q30.l125.fastq.gz \
    --outdir $CP/PLOIDY/$species \
    --mapper minimap2 \
    --preset sr \
    --caller varscan2 \
    --type both \
    --var $varscan \
    --min-var-freq 0.1 \
    --threads 12 \
    --mem 8
done

## Masked genomes
for species in {RCC138,RCC2335,RCC1871,CCMP1205,CCMP1998}; do 
  get_SNPs.pl \
    --fasta $CP/FASTA/Masked/$species.masked.fasta \
    --pe1 $CP/READS/FASTP/${species}_R1.q30.l125.fastq.gz \
    --pe2 $CP/READS/FASTP/${species}_R2.q30.l125.fastq.gz \
    --outdir $CP/PLOIDY/MASKED/$species \
    --mapper minimap2 \
    --preset sr \
    --caller varscan2 \
    --type both \
    --var $varscan \
    --min-var-freq 0.1 \
    --threads 12 \
    --mem 8
done

## not keeping the bam files with --bam, files are quite large

##### Plotting the ploidy as inferred from the allelic frequencies
## Unmasked
for species in {RCC138,RCC2335,RCC1871,CCMP1205,CCMP1998}; do 
  $CP/Scripts/split_VCF.pl \
    -v $CP/PLOIDY/$species/minimap2.varscan2.VCFs/*.vcf \
    -o $CP/PLOIDY/$species/VCFs

  $CP/Scripts/sort_SNPs.pl \
    -min 10 \
    -max 90 \
    -vcf $CP/PLOIDY/$species/VCFs/*.vcf \
    -noindel \
    -out $CP/PLOIDY/$species/TSV

  $CP/Scripts/plot_r.pl \
    -t $CP/PLOIDY/$species/TSV/*.tsv \
    -o $CP/PLOIDY/$species/R_plots \
    -ymax 0.12
done

## Masked (check effect on background noise for aneuploid chromosomes...)
# quick check for masked/unmasked metrics 
for species in {RCC138,RCC2335,RCC1871,CCMP1205,CCMP1998}; do
  $CP/Scripts/sort_SNPs.pl \
    -min 10 \
    -max 90 \
    -vcf $CP/PLOIDY/Masked/$species/minimap2.varscan2.VCFs/*.vcf \
    -noindel \
    -out $CP/PLOIDY/Masked/$species/TSV
done


## R plots (masked)
for species in {RCC138,RCC2335,RCC1871,CCMP1205,CCMP1998}; do 
  $CP/Scripts/split_VCF.pl \
    -v $CP/PLOIDY/MASKED/$species/minimap2.varscan2.VCFs/*.vcf \
    -o $CP/PLOIDY/MASKED/$species/VCFs

  $CP/Scripts/sort_SNPs.pl \
    -min 10 \
    -max 90 \
    -vcf $CP/PLOIDY/MASKED/$species/VCFs/*.vcf \
    -noindel \
    -out $CP/PLOIDY/MASKED/$species/TSV

  $CP/Scripts/plot_r.pl \
    -t $CP/PLOIDY/MASKED/$species/TSV/*.tsv \
    -o $CP/PLOIDY/MASKED/$species/R_plots \
    -ymax 0.12
done


###############################################################
### Getting normalized sequencing depths (on unmasked genomes)
###############################################################

for species in {RCC138,RCC1871,CCMP1205}; do
  cd $CP/PLOIDY/$species/minimap2.varscan2.coverage
  $CP/Scripts/Coverage_to_Circos.pl \
    -f *.coverage \
    -o $species \
    -p 2 \
    -n \
    -s 5000 \
    -w 10000
  cp $species.cov $CP/PLOIDY/
done

## RCC2335, CCMP1998 are haploid
for species in {RCC2335,CCMP1998}; do
  cd $CP/PLOIDY/$species/minimap2.varscan2.coverage
  $CP/Scripts/Coverage_to_Circos.pl \
    -f *.coverage \
    -o $species \
    -p 1 \
    -n \
    -s 5000 \
    -w 10000
  cp $species.cov $CP/PLOIDY/
done

cat $CP/PLOIDY/{RCC138,CCMP1205,RCC2335,RCC1871,CCMP1998}.cov > $CP/PLOIDY/illumina.depth.cov

cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/CCMP1205_vs_RCC138/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/CCMP1205_vs_RCC2335/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/CCMP1205_vs_RCC1871/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/CCMP1205_vs_CCMP1998/CIRCOS/CONCATENATED/

cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC2335_vs_RCC138/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC2335_vs_CCMP1205/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC2335_vs_RCC1871/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC2335_vs_CCMP1998/CIRCOS/CONCATENATED/

cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC1871_vs_RCC138/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC1871_vs_CCMP1205/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC1871_vs_RCC2335/CIRCOS/CONCATENATED/
cp $CP/PLOIDY/illumina.depth.cov $CP/SYNY/RCC1871_vs_CCMP1998/CIRCOS/CONCATENATED/

### modified the circos configuration files manually to add seq. depth normalizations

circos \
  -conf $CP/SYNY/CCMP1205_vs_RCC2335/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile CCMP1205_vs_RCC2335.png

circos \
  -conf $CP/SYNY/CCMP1205_vs_RCC138/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile CCMP1205_vs_RCC138.png

circos \
  -conf $CP/SYNY/CCMP1205_vs_RCC1871/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile CCMP1205_vs_RCC1871.png

circos \
  -conf $CP/SYNY/CCMP1205_vs_CCMP1998/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile CCMP1205_vs_CCMP1998.png

circos \
  -conf $CP/SYNY/RCC2335_vs_CCMP1205/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC2335_vs_CCMP1205.png

circos \
  -conf $CP/SYNY/RCC2335_vs_CCMP1998/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC2335_vs_CCMP1998.png

circos \
  -conf $CP/SYNY/RCC2335_vs_RCC138/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC2335_vs_RCC138.png

circos \
  -conf $CP/SYNY/RCC2335_vs_RCC1871/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC2335_vs_RCC1871.png

circos \
  -conf $CP/SYNY/RCC1871_vs_CCMP1998/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC1871_vs_CCMP1998.png

circos \
  -conf $CP/SYNY/RCC1871_vs_CCMP1205/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC1871_vs_CCMP1205.png

circos \
  -conf $CP/SYNY/RCC1871_vs_RCC2335/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC1871_vs_RCC2335.png

circos \
  -conf $CP/SYNY/RCC1871_vs_RCC138/CIRCOS/CONCATENATED/concatenated.inverted.mod.conf \
  -outputdir $CP/SYNY/ \
  -outputfile RCC1871_vs_RCC138.png



####################################################################################
### Telomeres + assembly metrics (in .telomeres_summary.txt files)
####################################################################################

for species in {RCC138,RCC2335,CCMP1205,RCC1871,CCMP1998}; do
  $CP/Scripts/check_for_telomeres.pl \
    -f $CP/FASTA/$species.fasta \
    -o TELOMERES \
    -p TTTAGG \
    -m 5
done

## CCMP1205/RCC138 telomeric unit is a 6-mer TTTAGG
## RCC2335 unit is a 8-mer TTTTAGG
## Homopolymer issues? Basecalling issues?
## TTTAGG, TTTTAGG and TTTTTAGG are found in all genomes

## Getting metrics from the GBFF files
$CP/Scripts/gbff_parser.pl \
  -i $CP/GBFF/*.gbff \
  -o $CP/METRICS


####################################################################################
### OrthoFinder 2.5.5
####################################################################################

orthofinder.py \
  -t 12 \
  -f $CP/OrthoFinder/FASTA \
  -o Results

## Creating Venn diagram from results

OGRES=/mnt/d/LinuxVM/Chloropicon/OrthoFinder/Results/Results_Nov16/Orthogroups

## Not sure why but some TSV files generated by OrthoFinder have \r\n line ends
## Even when generated on Linux; converting them with dos2unix

dos2unix $OGRES/*

$CP/Scripts/parse_OGs.pl \
  -og $OGRES/Orthogroups.tsv \
  $OGRES/Orthogroups_UnassignedGenes.tsv \
  -out $CP/OrthoFinder/Parsed_lists

## Getting unique orthogroups per species (excluding the other strains)
$CP/Scripts/get_species_set.pl \
  -l $CP/OrthoFinder/Parsed_lists/*.txt \
  -f $CP/OrthoFinder/Results/Results_Nov16/Orthogroup_Sequences \
  -d $CP/OrthoFinder/Results/Results_Nov16/Orthogroups/Orthogroups.tsv \
     $CP/OrthoFinder/Results/Results_Nov16/Orthogroups/Orthogroups_UnassignedGenes.tsv \
  -a $CP/SYNY/ALL/ANNOTATIONS/*.annotations \
  -x RCC138 CCMP1205 \
  -o $CP/OrthoFinder/C_roscoffensis \
  -p C_roscoffensis

$CP/Scripts/get_species_set.pl \
  -l $CP/OrthoFinder/Parsed_lists/*.txt \
  -f $CP/OrthoFinder/Results/Results_Nov16/Orthogroup_Sequences \
  -d $CP/OrthoFinder/Results/Results_Nov16/Orthogroups/Orthogroups.tsv \
     $CP/OrthoFinder/Results/Results_Nov16/Orthogroups/Orthogroups_UnassignedGenes.tsv \
  -a $CP/SYNY/ALL/ANNOTATIONS/*.annotations \
  -x RCC1871 RCC2335 CCMP1998 \
  -o $CP/OrthoFinder/C_primus \
  -p C_primus

