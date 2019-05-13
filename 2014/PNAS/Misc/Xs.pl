#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Calculates the number of unceartain amino acids (X), if present in protein fasta files

use strict;
use warnings;

my $usage = 'perl Xs.pl *.fasta';
die $usage unless @ARGV;

open OUT, ">protein_with_Xs.txt";

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fasta$//;
	$file =~ s/.fsa$//;
	my $protein = undef;
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>/ ){next;}
		else{$protein .= $line;}
	}
	if ($protein =~ /X/){
		print OUT "$file\n";
	}
}