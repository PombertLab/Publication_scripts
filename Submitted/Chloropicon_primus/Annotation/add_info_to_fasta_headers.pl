#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = add_info_to_fasta_headers.pl *.fasta"; die "\n$usage\n\n" unless @ARGV;

## NCBI Fasta headers 
my $organism = 'Prasinophyceae sp. CCMP1205';
my $strain = 'CCMP1205';
my $lineage = 'cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;';
my $gcode = '1';
my $moltype = 'genomic';
my $chromosome = 0;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.headers";
	$chromosome++;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			print OUT ">$1 [organism=$organism][strain=$strain][lineage=$lineage][gcode=$gcode][moltype=$moltype][chromosome=$chromosome]\n";
		}
		else {print OUT "$line\n";}
	}
	close OUT;
	system "mv $file.headers $file"; ## Overwrites original file 
}