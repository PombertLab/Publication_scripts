#!/usr/bin/perl

use strict;
use warnings;

my $usage = "perl parse_Uniques.pl min_evalue *.blastp *.tblastn";
die "$usage\n" unless @ARGV;

my $evalue = shift@ARGV;
my $scifmt = sprintf("%.100f", $evalue);

my %hash = ();
## Filling the hash
while (my $file = shift@ARGV){
	open IN, "<$file";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/){
			my $locus = $1; my $hit = $2; my $eval = $3;
			my $ev = sprintf("%.100f", $eval);
			if ($ev <= $scifmt){
				$hash{$locus} .= $hit;
			}
		}
	}
}
## Looking for uniques
open IN2, "<list.txt";
open OUT, ">not_unique.txt";
open OUT2, ">uniques.txt";
while (my $line = <IN2>){
	chomp $line;
	if (exists $hash{$line}){
		print OUT "$line\n";
	}
	else{
		print OUT2 "$line\n";
	}
}