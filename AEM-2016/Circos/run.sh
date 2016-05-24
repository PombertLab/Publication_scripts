#!/usr/bin/bash

perl 1_splitATCC_VCF.pl ATCC.vcf
perl 2_makeATCCbins.pl
perl 3_SNPs_slide.pl *.bin
perl 4_SNPsToCircos.pl 1000 c *.windows
cat *.SNPs > 1961-4_1000.SNPs