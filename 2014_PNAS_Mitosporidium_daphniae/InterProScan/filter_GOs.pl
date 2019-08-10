#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Parses the output of filter_InterProScan.pl
## Generates two lists: One to be used as input for cateGOrizer, the other one a simple list of the proteins and their GOs

use strict;
use warnings;

my $usage = 'perl filter_GOs.pl *.interpro.GOs'; die $usage unless@ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.out"; ## GOs per protein, duplicates filtered out
	open TEST, ">$file.categorizer"; ## Input list for CateGOrizer (http://www.animalgenome.org/bioinfo/tools/catego/) 
	while (my $line = <IN>){
		chomp $line;
		my %hash = ();
		if ($line =~ /^(DI\d+_\d+p\d+)\t(.*$)/){
			print OUT "$1\t";
			my @array = split(' ', $2);
			while (my $go = shift@array){
				if (exists $hash{$go}){next;}
				else {print OUT "$go\t"; print TEST "$go\n";$hash{$go}=$go;}
			}
			print OUT "\n";
		}
	}
}