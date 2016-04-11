#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl fasta_to_string.pl *.fsa';

die $usage unless @ARGV;

## convert the content of single fasta files to single & one-line DNA strings ##

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
