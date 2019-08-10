#!/usr/bin/perl

## convert the content of single fasta files to single & one-line DNA strings ##

use strict;
use warnings;

my $usage = 'perl fasta_to_string.pl *.fsa';

die $usage unless @ARGV;

while (my $fsa = shift @ARGV) {

	open IN, "<$fsa" or die "cannot open $fsa";
	$fsa =~ s/\.fsa$//;
	open OUT, ">$fsa.string";
	
	my @string = ();

	while (my $line = <IN>) {

		chomp $line;

		if ($line =~ /^>.*$/) {
				next;
		}
		else {
			push (@string, $line);
		}
		
		my $DNA = join('', @string);
		$DNA =~ s/\s//; 
		print OUT "$DNA";
		@string = ();
	}
	
	close IN;
}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)