#!/usr/bin/bash

## Grabbing read metrics with keep_longest_reads.pl (KLR) for MultiQC

## illumina data
~/GitHub/Misc/keep_longest_reads.pl \
  -o KLR \
  -p KLR_short_metrics \
  -x \
  -t \
  -j \
  -h 'Short read data' \
  -q 'KLR_SR' \
  -i E_cuniculi_50602/E_cuniculi_50602_S* \
  -i E_hellem_50604/E_hellem_50604_R* \
  -i E_intestinalis_50506/E_intestinalis_50506_R*

## nanopore + pacbio
~/GitHub/Misc/keep_longest_reads.pl \
  -o KLR \
  -p KLR_long_metrics \
  -x \
  -m 10000 \
  -t \
  -j \
  -h 'Long read data' \
  -q 'KLR_LR' \
  -i E_cuniculi_50602/E_cuniculi_50602.nanopore.fastq.gz \
  -i E_hellem_50604/E_hellem_50604.nanopore.fastq.gz \
  -i E_intestinalis_50506/E_intestinalis_50506.nanopore.fastq.gz \
  -i E_intestinalis_50506/E_intestinalis_50506_pacbio.fastq.gz

