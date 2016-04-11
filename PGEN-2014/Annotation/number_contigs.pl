#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl number_contigs.pl *.fsa';

die $usage unless @ARGV;

## Change the name of the contigs so that they increment by 1 ##

my $counter = 1;

while (my $fsa = shift @ARGV) {

	open IN, "<$fsa" or die "cannot open $fsa";
	$fsa =~ s/(d+)\.fsa$//;
	open OUT, ">contig-$counter.fasta";
	
	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^>contig-.*$/) {
			print OUT ">contig-$counter.fasta\n";
		}
			
		else {
			print OUT "$line\n";
		}
	}
	
	$counter++;
	
	close IN;
}

exit;
