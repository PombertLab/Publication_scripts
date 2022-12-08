#!/usr/bin/bash

# Example command to output basecalls, mappings, and CpG 5mC and 5hmC methylation in both per-read (``mod_mappings``)
# and aggregated (``mods``) formats
# Compute settings: GPU devices 0 with 12 CPU cores
# IO intensive, trying from SSD
# 

## Run in conda: conda activate megalodon

# Data location
FAST5=/media/FatCat/Nanopore/E_cuniculi_ATCC_50602/Cuni/20201124_1952_MN23338_FAO08622_6cc86583/fast5/
FASTA=/media/Plab/E_cuniculi_50602/Cuniculi_50602_v3.fasta
OUTDIR=/media/Plab/E_cuniculi_50602/Megalodon_hac/

# Model parameters
# To list possible remora models:
# remora model list_pretrained
BCMODEL=dna_r9.4.1_450bps_hac.cfg
RMMODEL=dna_r9.4.1_e8
RMTYPE=hac
GUPPYSERV=/media/Data_1/opt/ont-guppy/bin/guppy_basecall_server

## Make sure to update remora model to hac if hac and fast if fast
## Change model to e8.1 (latest); Can't e8.1 hac only detects 5mc not 5hmc
## Must stick to e8.0
megalodon \
    $FAST5 \
    --guppy-config $BCMODEL \
    --remora-modified-bases $RMMODEL $RMTYPE 0.0.0 5hmc_5mc CG 0 \
    --output-directory $OUTDIR \
    --outputs basecalls mappings mod_mappings mods \
    --write-mod-log-probs \
	--sort-mappings \
    --mod-min-prob 0.1 \
    --reference $FASTA \
    --devices 0 \
    --processes 12 \
    --guppy-server-path $GUPPYSERV

## --sort-mappings sorts + index the bam file

