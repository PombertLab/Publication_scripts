#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'run_phyloFit_ORDOSPORA';
my $version = 0.1;

use strict; use warnings;

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Runs phyloFit and phastBias from the PHAST package (http://compgen.cshl.edu/phast/) on datasets prepared with make_PHAST_datasets.pl	
EXAMPLE		run_phyloFit_ORDOSPORA.pl CWI40 *.NT.macse
NOTE		This script is hard-coded for the Ordospora dataset 
OPTIONS
die "$options\n" unless @ARGV;

my $branch = shift@ARGV;
while (my $file = shift @ARGV){
	print "$file\n";
	system "phyloFit --tree \"((CWI40,CWI41),(CWI42,M896))\" -o $file.$branch $file";
	system "phastBias $file $file.$branch.mod $branch --posteriors wig --output-tracts $file.$branch.gBGC.gff > $file.$branch.scores.wig";
}