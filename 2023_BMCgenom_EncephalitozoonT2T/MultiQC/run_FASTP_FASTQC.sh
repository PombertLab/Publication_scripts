#!/usr/bin/bash

# ln -s ../illumina/C1_S1_L001_R1_001.fastq.gz E_cuniculi_50602_S1_R1.fastq.gz
# ln -s ../illumina/C1_S1_L001_R2_001.fastq.gz E_cuniculi_50602_S1_R2.fastq.gz
# ln -s ../illumina/C2_S2_L001_R1_001.fastq.gz E_cuniculi_50602_S2_R1.fastq.gz
# ln -s ../illumina/C2_S2_L001_R2_001.fastq.gz E_cuniculi_50602_S2_R2.fastq.gz

S1_R1=E_cuniculi_50602_S1_R1.fastq.gz
S1_R2=E_cuniculi_50602_S1_R2.fastq.gz
S2_R1=E_cuniculi_50602_S2_R1.fastq.gz
S2_R2=E_cuniculi_50602_S2_R2.fastq.gz

fastp \
  -w 10 \
  -i $S1_R1 \
  -I $S1_R2 \
  -o E_cuniculi_50602_S1_R1.q30.l125.fastq.gz \
  -O E_cuniculi_50602_S1_R2.q30.l125.fastq.gz \
  -M 30 \
  -r \
  -l 125 \
  --adapter_sequence=CTGTCTCTTATACACATCT \
  --adapter_sequence_r2=CTGTCTCTTATACACATCT \
  -j E_cuniculi_50602_S1_fastp.json \
  -h E_cuniculi_50602_S1_fastp.html 

fastp \
  -w 10 \
  -i $S2_R1 \
  -I $S2_R2 \
  -o E_cuniculi_50602_S2_R1.q30.l125.fastq.gz \
  -O E_cuniculi_50602_S2_R2.q30.l125.fastq.gz \
  -M 30 \
  -r \
  -l 125 \
  --adapter_sequence=CTGTCTCTTATACACATCT \
  --adapter_sequence_r2=CTGTCTCTTATACACATCT \
  -j E_cuniculi_50602_S2_fastp.json \
  -h E_cuniculi_50602_S2_fastp.html 

fastqc *.gz
