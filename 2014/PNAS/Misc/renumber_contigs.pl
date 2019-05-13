#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Renames contigs so that they increment by 1

use strict;
use warnings;

my $usage = 'perl number_contigs.pl *.fasta';
die $usage unless @ARGV;

my $counter = 1;

while (my $fasta = shift@ARGV) {

	open IN, "<$fasta" or die "cannot open $fasta";
	$fasta =~ s/.fasta$//;
	$fasta =~ s/.fsa$//;
	open OUT, ">$fasta.input";
	
	while (my $line = <IN>) {
		chomp $line;
		if ($line =~ /^>/) {
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
