#!/usr/bin/env bash

GENDIV=/home/jpombert/Analyses/Encephalitozoon/GenDiv
mkdir -p $GENDIV/FASTQ/FASTQ_ILLUMINA
cd $GENDIV/FASTQ/FASTQ_ILLUMINA

## NCBI SRA toolkit version => sratoolkit.3.0.1

##########   Encephalitozoon hellem   ##########

fasterq-dump SRR17853475 --split-3 ## Encephalitozoon hellem ATCC 50604
fasterq-dump SRR23560257 --split-3 ## Encephalitozoon hellem ATCC 50451 set 1
fasterq-dump SRR23560258 --split-3 ## Encephalitozoon hellem ATCC 50451 set 2
fasterq-dump SRR14062087 --split-3 ## Encephalitozoon hellem Swiss Isolate set1, 54 bp GAIIx
fasterq-dump SRR14017862 --split-3 ## Encephalitozoon hellem Swiss Isolate set2, 54 bp GAIIx

for x in {SRR17853475,SRR14062087,SRR23560257,SRR23560258,SRR14017862}; do
	mv ${x}_1.fastq ${x}.R1.fastq;
	mv ${x}_2.fastq ${x}.R2.fastq;
done

mv SRR17853475.R1.fastq 50604.raw.R1.fastq
mv SRR17853475.R2.fastq 50604.raw.R2.fastq

cat SRR23560257.R1.fastq SRR23560258.R1.fastq > 50451.raw.R1.fastq
cat SRR23560257.R2.fastq SRR23560258.R2.fastq > 50451.raw.R2.fastq
rm SRR23560257* SRR23560258*

cat SRR14062087.R1.fastq SRR14017862.R1.fastq > Swiss.raw.R1.fastq
cat SRR14062087.R2.fastq SRR14017862.R2.fastq > Swiss.raw.R2.fastq
rm SRR14062087* SRR14017862*

##########   Encephalitozoon intestinalis   ##########

fasterq-dump SRR17865591 --split-3 ## Encephalitozoon intestinalis ATCC 50506
fasterq-dump SRR24007516 --split-3 ## Encephalitozoon intestinalis ATCC 50507
fasterq-dump SRR24007515 --split-3 ## Encephalitozoon intestinalis ATCC 50603
fasterq-dump SRR24007514 --split-3 ## Encephalitozoon intestinalis ATCC 50651

for x in {SRR17865591,SRR24007514,SRR24007515,SRR24007516}; do
	mv ${x}_1.fastq ${x}.R1.fastq;
	mv ${x}_2.fastq ${x}.R2.fastq;
done

mv SRR17865591_1.fastq 50506.raw.R1.fastq
mv SRR17865591_2.fastq 50506.raw.R2.fastq

mv SRR24007516_1.fastq 50507.raw.R1.fastq
mv SRR24007516_2.fastq 50507.raw.R2.fastq

mv SRR24007515_1.fastq 50603.raw.R1.fastq
mv SRR24007515_2.fastq 50603.raw.R2.fastq

mv SRR24007514_1.fastq 50651.raw.R1.fastq
mv SRR24007514_2.fastq 50651.raw.R2.fastq

##########   Encephalitozoon cuniculi   ##########

## Encephalitozoon cuniculi ATCC 50602
fasterq-dump SRR17858635 --split-3 
fasterq-dump SRR17858636 --split-3

cat SRR17858635_1.fastq SRR17858636_1.fastq > 50602.raw.R1.fastq
cat SRR17858635_2.fastq SRR17858636_2.fastq > 50602.raw.R2.fastq
pigz 50602.raw.R1.fastq
pigz 50602.raw.R2.fastq

rm SRR17858635_1.fastq SRR17858635_2.fastq SRR17858636_1.fastq SRR17858636_2.fastq

## Encephalitozoon cuniculi ECI
fasterq-dump SRR122312 --split-3 ## Encephalitozoon cuniculi ECI
fasterq-dump SRR122313 --split-3 ## Encephalitozoon cuniculi ECI
fasterq-dump SRR122314 --split-3 ## Encephalitozoon cuniculi ECI

cat SRR122312_1.fastq SRR122313_1.fastq SRR122314_1.fastq > ECI.raw.R1.fastq
cat SRR122312_2.fastq SRR122313_2.fastq SRR122314_2.fastq > ECI.raw.R2.fastq
pigz ECI.raw.R1.fastq
pigz ECI.raw.R2.fastq

rm SRR122312_1.fastq SRR122313_1.fastq SRR122314_1.fastq
rm SRR122312_2.fastq SRR122313_2.fastq SRR122314_2.fastq

## Encephalitozoon cuniculi ECII
fasterq-dump SRR122309 --split-3
fasterq-dump SRR122311 --split-3
fasterq-dump SRR122315 --split-3

cat SRR122309_1.fastq SRR122311_1.fastq SRR122315_1.fastq > ECII.raw.R1.fastq
cat SRR122309_2.fastq SRR122311_2.fastq SRR122315_2.fastq > ECII.raw.R2.fastq
pigz ECII.raw.R1.fastq
pigz ECII.raw.R2.fastq

rm SRR122309_1.fastq SRR122311_1.fastq SRR122315_1.fastq
rm SRR122309_2.fastq SRR122311_2.fastq SRR122315_2.fastq

## Encephalitozoon cuniculi ECII-CZ; single not paired ends!
fasterq-dump SRR636827 --split-3
mv SRR636827.fastq ECII_CZ.raw.SE.fastq
pigz ECII_CZ.raw.SE.fastq

## Encephalitozoon cuniculi ECIII
fasterq-dump SRR122310 --split-3
fasterq-dump SRR122316 --split-3
fasterq-dump SRR122317 --split-3

cat SRR122310_1.fastq SRR122316_1.fastq SRR122317_1.fastq > ECIII.raw.R1.fastq
cat SRR122310_2.fastq SRR122316_2.fastq SRR122317_2.fastq > ECIII.raw.R2.fastq

pigz ECIII.raw.R1.fastq
pigz ECIII.raw.R2.fastq

rm SRR122310_1.fastq SRR122316_1.fastq SRR122317_1.fastq
rm SRR122310_2.fastq SRR122316_2.fastq SRR122317_2.fastq

## Encephalitozoon cuniculi ECIII-L
fasterq-dump SRR2105612 --split-3
mv SRR2105612_1.fastq ECIII_L.raw.R1.fastq
mv SRR2105612_2.fastq ECIII_L.raw.R2.fastq

pigz ECIII_L.raw.R1.fastq
pigz ECIII_L.raw.R2.fastq


##########   Ordospora colligata   ##########

# NO-V-7
fasterq-dump SRR9597068	--split-3
mv SRR9597068_1.fastq NOV7.raw.R1.fastq
mv SRR9597068_2.fastq NOV7.raw.R2.fastq
pigz NOV7.raw.R1.fastq
pigz NOV7.raw.R2.fastq

# GB-EP-1
fasterq-dump SRR9597069	--split-3
mv SRR9597069_1.fastq GBEP.raw.R1.fastq
mv SRR9597069_2.fastq GBEP.raw.R2.fastq
pigz GBEP.raw.R1.fastq
pigz GBEP.raw.R2.fastq

# FI-SK-17-1
fasterq-dump SRR9597070	--split-3
mv SRR9597070_1.fastq FISK.raw.R1.fastq
mv SRR9597070_2.fastq FISK.raw.R2.fastq
pigz FISK.raw.R1.fastq
pigz FISK.raw.R2.fastq

