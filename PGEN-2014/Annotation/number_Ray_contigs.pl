#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl number_contigs.pl *.fasta';

die $usage unless @ARGV;

## Change the name of the contigs produced by Ray so that they increment by 1 ##

my $counter = 1;

while (my $fasta = shift @ARGV) {

	open IN, "<$fasta" or die "cannot open $fasta";
	$fasta =~ s/.fasta$//;
	open OUT, ">$fasta.fsa";
	
	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^>contig-.*$/) {
			print OUT ">contig-$counter\n";
			$counter++;
		}
			
		else {
			print OUT "$line\n";
		}
	}
	
	close IN;
}

exit;
