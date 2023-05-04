#!/usr/bin/bash

#### https://github.com/PacificBiosciences/pbbioconda

conda create --name pacbio
conda activate pacbio

conda install -c bioconda pbccs
conda install -c bioconda primrose
conda install -c bioconda pbcommand
conda install -c bioconda pbmm2
conda install -c bioconda pbbam

##### Run CCS inside conda pacbio environment
ROOT=/media/Data_2/Sequences/E_intestinalis_ATCC50506_PacBio
SUB=isilon/seq/smrt3/data_root/runs/54027/r54027_20190906_191401/1_A01
BAM=m54027_190906_192305.subreads.bam

ccs \
  $ROOT/$SUB/$BAM \
  $ROOT/E_intestinalis.hifi_reads.bam \
  --hifi-kinetics

##### Run primrose inside conda pacbio environment
primrose \
  $ROOT/E_intestinalis.hifi_reads.bam \
  $ROOT/E_intestinalis.5mC.hifi_reads.bam


## Splitting Multifasta to single fasta
cd /media/Data_2/Sequences/E_intestinalis_ATCC50506_PacBio/Primrose/

split_Fasta.pl \
  -f GCA_024399295.1_ASM2439929v1_genomic.fasta \
  -o SPLIT_FASTA

cd /media/Data_2/Sequences/E_intestinalis_ATCC50506_PacBio/Primrose/SPLIT_FASTA

## Primrose v1.3.0 - https://github.com/PacificBiosciences/primrose
# Aligning to 1 chromosome at a time to prevent issue where all reads
# from rRNA map to the first instance encountered when mapped with pbmm2

for ACC in {CP075158,CP075159,CP075160,CP075161,CP075162,CP075163,CP075164,CP075165,CP075166,CP075167,CP075168}; do 
  pbmm2 index $ACC.fasta $ACC.mmi; 
  pbmm2 align $ACC.mmi ../E_intestinalis.5mC.hifi_reads.bam $ACC.ref.bam --sort -j 4 -J 2;
  pbindex $ACC.ref.bam;
done

conda deactivate

## CpG tools v1.1.0 - https://github.com/PacificBiosciences/pb-CpG-tools
# git clone https://github.com/PacificBiosciences/pb-CpG-tools.git
# conda env create -f conda_env_cpg.yaml
conda activate cpg

## Note: struggles if bam file is in a different location...
## looks sensitive to pathing issues

for ACC in {CP075158,CP075159,CP075160,CP075161,CP075162,CP075163,CP075164,CP075165,CP075166,CP075167,CP075168}; do 
  ../pb-CpG-tools/aligned_bam_to_cpg_scores.py \
    -b $ACC.ref.bam \
    -m reference \
    -f $ACC.fasta \
    -o $ACC \
    -d ../pb-CpG-tools/pileup_calling_model/;
done
