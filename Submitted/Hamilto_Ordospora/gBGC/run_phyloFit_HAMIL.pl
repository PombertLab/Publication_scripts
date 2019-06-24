#!/usr/bin/perl

use strict;
use warnings;

die "\nUSAGE = run_phyloFit.pl CWI36 *.NT.macse\n\n" unless @ARGV;

my $branch = shift@ARGV;
while (my $file = shift @ARGV){
	print "$file\n";
	system "phyloFit --tree \"((CWI36,CWI39),(CWI37,CWI38))\" -o $file.$branch $file";
	system "phastBias $file $file.$branch.mod $branch --posteriors wig --output-tracts $file.$branch.gBGC.gff > $file.$branch.scores.wig";
}