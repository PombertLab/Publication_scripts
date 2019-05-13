#!/usr/bin/perl

## Catenate all SNPs metrics from the GeneAlignment_SNPs.pl script into one file

use strict;
use warnings;

my $usage = 'perl cat_metrics.pl *.metrics';

die $usage unless @ARGV;

open OUT, ">all_gene_SNPs_metrics.txt"; ## Change to desired file name
print OUT "Gene\tlength\tSNPs\tRatio\n";

while (my $metrics = shift @ARGV) {

	open IN, "<$metrics" or die "cannot open $metrics";
	
		while (my $line = <IN>) {

			chomp $line;
		
			if ($line =~ /^ECU/) {
			
			print OUT "$line\n";
			
			}
		}
	close IN;
}

close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)