#!/usr/bin/bash

# Fedora Linux 35 (Workstation Edition
# AMD Ryzen 9 3900X (24) @ 3.800GHz
# NVIDIA RTX A4000

## Install tombo (worked only with python 3.7)
## conda create --name tombo python=3.7
## conda activate tombo
## conda install -c bioconda ont-tombo
## conda install numpy=1.19 ## solves ValueError: cannot convert float NaN to integer with numpy >= 1.20

###### Tombo DOES NOT WORK on multi fast5 files !!!!! #####
## Requires basecalled, single fast5:
## 1. Basecalled everything again with guppy
## 2. Split the fast5 files

################## Basecall with guppy #####################
FAST5=/media/FatCat/Nanopore/E_cuniculi_ATCC_50602/Cuni/20201124_1952_MN23338_FAO08622_6cc86583/fast5
OUTPUT=/media/Plab/E_cuniculi_50602/Tombo/basecall
LOG=/media/Plab/E_cuniculi_50602/Tombo/basecall/guppy.log

## basecall with fast5
echo "Started on `date`" >> $LOG
guppy_basecaller \
	-i $FAST5 \
	-r \
	-s $OUTPUT \
	--fast5_out \
	--bam_out \
	--flowcell FLO-MIN106 \
	--kit SQK-LSK109 \
	--trim_adapters \
	-x 'cuda:0' 
echo "Completed on `date`" >> $LOG

##### Split FAST5 files with multi_to_single_fast5
## pip install ont-fast5-api or use the one from Megalodon
F5=/media/Plab/E_cuniculi_50602/Tombo/basecall/workspace
SINGLES=/media/Plab/E_cuniculi_50602/Tombo/f5_single

/home/jpombert/.conda/envs/megalodon/bin/multi_to_single_fast5 \
  -i $F5 \
  -s $SINGLES \
  -t 12

## Create read subset; use mv to save SSD space
SUBSET=/media/Plab/E_cuniculi_50602/Tombo/reads_subset/
mkdir -p $SUBSET
for k in {0..124}; do mv $SINGLES/$k $SUBSET; done


###################### Tombo #########################
FAST5=$SUBSET
FASTA=/media/Plab/E_cuniculi_50602/Cuniculi_50602_v3.fasta
NAME=ecun_50602

## Requires annotated fast5 files; --overwrites replaces the tombo index files if present
tombo resquiggle \
  $FAST5 \
  $FASTA \
  --processes 12 \
  --num-most-common-errors 5 \
  --overwrite

############## De novo
## de novo; # works with 40k reads, crashes with 500K,
## hard drives trashes, move to SSD? -> Works on SSD
tombo detect_modifications de_novo \
  --fast5-basedirs $FAST5 \
  --statistics-file-basename $NAME.de_novo_detect \
  --processes 12

tombo text_output browser_files \
  --statistics-filename $NAME.de_novo_detect.tombo.stats \
  --browser-file-basename $NAME.de_novo_detect \
  --file-types dampened_fraction


## https://nanoporetech.github.io/tombo/text_output.html
# de novo
tombo text_output signif_sequence_context \
  --statistics-filename $NAME.de_novo_detect.tombo.stats \
  --genome-fasta $FASTA \
  --num-regions 1000 \
  --num-bases 50

/media/Data_1/opt/meme/bin/meme \
  -oc tombo.de_novo_motif_detection.meme \
  -dna \
  -mod zoops \
  tombo_results.significant_regions.fasta

# plot raw signal at most significant denovo locations
tombo plot most_significant \
  --fast5-basedirs $FAST5 \
  --statistics-filename $NAME.de_novo_detect.tombo.stats \
  --plot-standard-model \
  --pdf-filename $NAME.most_significant_denovo_sites.pdf

MEMEDIR=/media/Plab/E_cuniculi_50602/Tombo/tombo.de_novo_motif_detection.meme/
mkdir -p $MEMEDIR
mv tombo_results.significant_regions.fasta $MEMEDIR

############## dam/dcm
tombo detect_modifications alternative_model \
  --fast5-basedirs $FAST5 \
  --statistics-file-basename $NAME \
  --alternate-bases dam dcm \
  --processes 12

# produces wig file with estimated fraction of modified reads at each valid reference site
tombo text_output browser_files \
  --statistics-filename $NAME.dam.tombo.stats \
  --file-types dampened_fraction \
  --browser-file-basename $NAME.dam

tombo text_output browser_files \
  --statistics-filename $NAME.dcm.tombo.stats \
  --file-types dampened_fraction \
  --browser-file-basename $NAME.dcm

# also produce successfully processed reads coverage file for reference
tombo text_output browser_files \
  --fast5-basedirs $FAST5 \
  --file-types coverage \
  --browser-file-basename $NAME

# plot raw signal at most significant dam locations
tombo plot most_significant \
  --fast5-basedirs $FAST5 \
  --statistics-filename $NAME.dam.tombo.stats \
  --plot-standard-model \
  --plot-alternate-model dam \
  --pdf-filename $NAME.most_significant_dam_sites.pdf

# plot raw signal at most significant dcm locations
tombo plot most_significant \
  --fast5-basedirs $FAST5 \
  --statistics-filename $NAME.dcm.tombo.stats \
  --plot-standard-model \
  --plot-alternate-model dcm \
  --pdf-filename $NAME.most_significant_dcm_sites.pdf

# dam
tombo text_output signif_sequence_context \
  --statistics-filename $NAME.dam.tombo.stats \
  --genome-fasta $FASTA \
  --num-regions 1000 \
  --num-bases 50

/media/Data_1/opt/meme/bin/meme \
  -oc tombo.dam_detection.meme \
  -dna \
  -mod zoops \
  tombo_results.significant_regions.fasta

MEMEDIR=/media/Plab/E_cuniculi_50602/Tombo/tombo.dam_detection.meme/
mkdir -p $MEMEDIR
mv tombo_results.significant_regions.fasta $MEMEDIR

# dcm
tombo text_output signif_sequence_context \
  --statistics-filename $NAME.dcm.tombo.stats \
  --genome-fasta $FASTA \
  --num-regions 1000 \
  --num-bases 50

/media/Data_1/opt/meme/bin/meme \
  -oc tombo.dcm_detection.meme \
  -dna \
  -mod zoops \
  tombo_results.significant_regions.fasta

MEMEDIR=/media/Plab/E_cuniculi_50602/Tombo/tombo.dcm_detection.meme/
mkdir -p $MEMEDIR
mv tombo_results.significant_regions.fasta $MEMEDIR
