#!/usr/bin/bash

# Example command to output basecalls, mappings, and CpG 5mC and 5hmC methylation in both per-read (``mod_mappings``)
# and aggregated (``mods``) formats
# Compute settings: GPU devices 0 with 12 CPU cores

## pip install megalodon==2.5.0
## wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.4.6_linux64.tar.gz

## Looks like the 5hmC model disappeared in the new remora, run with older version??
## pip install megalodon==2.5.0
## pip install ont-remora==2.0.0
## pip index versions ont-remora

## On Milady
## Guppy

# Data location
FAST5=/media/FatCat_2/Microsporidia/Original_Nanopore_Data/E_hellem_ATCC50451/HEL50451/20201126_1841_MN23338_FAO94202_7a90a7e0/fast5
FASTA=/home/jpombert/FASTA/Hellem_T2T_ATCC_50451.fna
OUTDIR=/media/FatCat_1/jpombert/Megalodon/50451/

# Model parameters
# To list possible remora models:
# remora model list_pretrained
BCMODEL=dna_r9.4.1_450bps_hac.cfg
RMMODEL=dna_r9.4.1_e8
RMTYPE=hac
GUPPY=/home/jpombert/ont-guppy-6.4.6//bin
GUPPYSERV=$GUPPY/guppy_basecall_server

## Must stick to e8.0; 5hmC no longer available in remora 2.0.0
megalodon \
    $FAST5 \
    --guppy-config $BCMODEL \
    --remora-modified-bases $RMMODEL $RMTYPE 0.0.0 5mc CG 0 \
    --output-directory $OUTDIR \
    --outputs basecalls mappings mod_mappings mods \
    --write-mod-log-probs \
    --sort-mappings \
    --mod-min-prob 0.1 \
    --reference $FASTA \
    --devices 0 \
    --processes 12 \
    --guppy-server-path $GUPPYSERV \
    --overwrite

## --sort-mappings sorts + index the bam file

