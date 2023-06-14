#!/usr/bin/env bash

GENDIV=/home/jpombert/Analyses/Encephalitozoon/GenDiv
R150=/home/jpombert/Analyses/Encephalitozoon/GenDiv/RM/SSRG_150/minimap2.varscan2.VCFs
R250=/home/jpombert/Analyses/Encephalitozoon/GenDiv/RM/SSRG_250/minimap2.varscan2.VCFs
ILLU=/home/jpombert/Analyses/Encephalitozoon/GenDiv/RM/ILLUMINA/minimap2.varscan2.VCFs

## E. hellem
# no illumina data for 50504
$GENDIV/Scripts/gceq.pl \
  -f $GENDIV/FASTA/{Swiss,50451,50504,50604}.fasta \
  -v \
  $R150/{50504,50451,Swiss,50604}.*.vcf \
  $R250/{50504,50451,Swiss,50604}.*.vcf \
  $ILLU/{50451,Swiss,50604}.*.vcf \
  -o $GENDIV/GCeq/E_hellem.gceq \
  -verbose

## E. intestinalis
$GENDIV/Scripts/gceq.pl \
  -f $GENDIV/FASTA/{50506,50507,50603,50651}.fasta \
  -v \
  $R150/{50506,50507,50603,50651}.*.vcf \
  $R250/{50506,50507,50603,50651}.*.vcf \
  $ILLU/{50506,50507,50603,50651}.*.vcf \
  -o $GENDIV/GCeq/E_intestinalis.gceq \
  -verbose

## E. cuniculi
# no illumina data for GBM1
$GENDIV/Scripts/gceq.pl \
  -f $GENDIV/FASTA/{50602,GB_M1,ECI,ECII,ECII_CZ,ECIII,ECIII_L}.fasta \
  -v \
  $R150/{50602,GB_M1,ECI,ECII,ECII_CZ,ECIII,ECIII_L}.*.vcf \
  $R250/{50602,GB_M1,ECI,ECII,ECII_CZ,ECIII,ECIII_L}.*.vcf \
  $ILLU/{50602,ECI,ECII,ECII_CZ,ECIII,ECIII_L}.*.vcf \
  -o $GENDIV/GCeq/E_cuniculi.gceq \
  -verbose

## O. colligata
# no illumina data for OC4
$GENDIV/Scripts/gceq.pl \
  -f $GENDIV/FASTA/{OC4,NOV7,FISK,GBEP}.fasta \
  -v \
  $R150/{OC4,NOV7,FISK,GBEP}.*.vcf \
  $R250/{OC4,NOV7,FISK,GBEP}.*.vcf \
  $ILLU/{NOV7,FISK,GBEP}.*.vcf \
  -o $GENDIV/GCeq/O_coligata.gceq \
  -verbose
