#!/usr/bin/perl

## Filters pairwise SNPs 
## Works with the .snps output of cat_gene_SNPs.pl

use strict;
use warnings;

my $usage = 'perl pairwise_gene_SNPs_I_II.pl *.snps';

die $usage unless @ARGV;

while (my $snps = shift @ARGV) {

	open IN, "<$snps" or die "cannot open $snps";
	$snps =~ s/\.snps$//;
	open OUT, ">$snps.I_II.snps_pairs";

	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^(ECU\d+_\S+)\t(\d+)\t(\w)\t(\w)\t(\w)/) {
		
		my $gene = $1;
		my $position = $2;
		my $strain1 = $3;
		my $strain2 = $4;
		my $strain3 = $5;

			if ($strain1 ne $strain2) {
			
			print OUT "$gene\t$position\t$strain1\t$strain2\n";
			
			}

		}
		
	}
	
	close IN;
	close OUT;
	
}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)