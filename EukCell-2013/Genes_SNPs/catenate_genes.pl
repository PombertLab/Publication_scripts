#!/usr/bin/perl

## Catenate single fasta files (ECI, ECII and ECIII) into mulitfasta files to be used as inputs for gene alignments
## The list provided should contains one gene per line

use strict;
use warnings;

my $usage = 'perl catenate_genes.pl *.txt';

die $usage unless @ARGV;

while (my $txt = shift @ARGV) {

	open IN, "<$txt" or die "cannot open $txt";
	
	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^(ECU\d+.*)/) {
		
		my $ECU = $1;
		
		system "cat ECI_$ECU.fsa ECII_$ECU.fsa ECIII_$ECU.fsa > $ECU.fasta";
		
		}
	}
	
	close IN;
}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)
