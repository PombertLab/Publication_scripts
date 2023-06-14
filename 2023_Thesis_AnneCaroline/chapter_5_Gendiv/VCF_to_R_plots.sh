#!/usr/bin/bash

## Creating chromosome core (CC) and subtelomere (SUB) VCF subsets
# nanopore
$GENDIV/Scripts/parse_VCF_by_CC.pl \
  -l $GENDIV/Scripts/cc_list.txt \
  -v $GENDIV/Nanopore/minimap2.varscan2.VCFs/{50451,50506,50602,50604}*.vcf \
  -o $GENDIV/Nanopore/CC

# illumina
$GENDIV/Scripts/parse_VCF_by_CC.pl \
  -l $GENDIV/Scripts/cc_list.txt \
  -v $GENDIV/RM/HZ/minimap2.varscan2.VCFs/{50451,50506,50602,50604}.*.vcf \
  -o $GENDIV/RM/HZ/CC

## Preparing distributions for plotting with R
$GENDIV/Scripts/sort_SNPs.pl \
  -noindel \
  -vcf \
  $GENDIV/Nanopore/CC/*.vcf \
  $GENDIV/Nanopore/minimap2.varscan2.VCFs/*.vcf \
  $GENDIV/RM/HZ/CC/*.vcf \
  $GENDIV/RM/HZ/minimap2.varscan2.VCFs/{50451,50506,50602,50604}.*.vcf \
  -o $GENDIV/for_R_plots

## Plotting with R
$GENDIV/Scripts/plot_r.pl \
  -t $GENDIV/for_R_plots/*.tsv \
  -o $GENDIV/for_R_plots/PDFs \
  -ymax 0.12
