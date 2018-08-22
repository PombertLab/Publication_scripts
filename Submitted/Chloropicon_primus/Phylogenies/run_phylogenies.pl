#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

die "\nUSAGE = run_phylogenies.pl -t 10 -f *.txt

-t (--threads)	Number of threads to use [Default: 10]
-f (--fasta)	Files in Multifasta format

## Requires, MAFFT, TRIMAL and IQTREE
\n" unless @ARGV;

my $threads = 10;
my @files;
GetOptions(
	't|threads=i'=>\$threads,
	'f|fasta=s@{1,}'=>\@files
);

while (my $fasta = shift@files){
	$fasta =~ s/.txt$//; $fasta =~ s/.fasta$//; $fasta =~ s/.fsa$//;
	system "echo Aligning $fasta with MAFTT...";
	system "mafft $fasta.txt > $fasta.aln";
	system "echo Trimming $fasta alignment with TrimAL...";
	system "trimal -in $fasta.aln -out ${fasta}_trimmed.tfa -htmlout $fasta.html -automated1";
	system "Computing $fasta phylogenetic tree with IQtree...";
	system "iqtree -s ${fasta}_trimmed.tfa -nt $threads -bb 1000";
}
