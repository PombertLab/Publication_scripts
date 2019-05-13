#!/usr/bin/perl

## Catenate all gene SNPs from the GeneAlignment_SNPs.pl script into one file

use strict;
use warnings;

my $usage = 'perl cat_gene_SNPs.pl *.snps';

die $usage unless @ARGV;

open OUT, ">all_gene_snps.snps"; ## Change to desired file name

while (my $snps = shift @ARGV) {

	open IN, "<$snps" or die "cannot open $snps";
	$snps =~ s/.snps//;
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^(\d+)\t(\w)\t(\w)\t(\w)/) {
			
			print OUT "$snps\t$line\n";
			
			}
		}
	close IN;
}

close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)