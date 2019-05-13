#!/usr/bin/perl

## Catenate all codon SNPs into one file

use strict;
use warnings;

my $usage = 'perl cat_codon_SNPs.pl *.snps';

die $usage unless @ARGV;

open OUT, ">all_codon_snps.snps";

while (my $snps = shift @ARGV) {

	open IN, "<$snps" or die "cannot open $snps";
	$snps =~ s/.codon.snps//;
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^(\d+)\t(\d+)\t(\w{3})\t(\w{3})\t(\w{3})/) {
			
			print OUT "$snps\t$line\n";
			
			}
		}
	close IN;
}

close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)