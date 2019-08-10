#!/usr/bin/perl

use warnings;
use strict;

## Rename Chlre3 protein catalogs with simpler fasta headers

my $usage = 'perl rename_chlamy_prots.pl proteins.fasta';
die $usage unless @ARGV;
my $fasta = $ARGV[0];

open IN, "<$fasta";
open OUT, ">$fasta.renamed";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>jgi\|Chlre3\|(\d{1,})/){
		my $gene = $1;
		print OUT ">$gene\n";
	}else{
		print OUT "$line\n";
	}
}

