#!/usr/bin/perl

use strict;
use warnings;

die "\nUSAGE = run_phyloFit.pl CWI36 *.NT.macse\n\n" unless @ARGV;

my $branch = shift@ARGV;
while (my $file = shift @ARGV){
	print "$file\n";
	system "phyloFit --tree \"((CWI40,CWI41),(CWI42,M896))\" -o $file.$branch $file";
	system "phastBias $file $file.$branch.mod $branch --posteriors wig --output-tracts $file.$branch.gBGC.gff > $file.$branch.scores.wig";
}