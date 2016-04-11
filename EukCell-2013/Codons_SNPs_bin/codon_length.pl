#!/usr/bin/perl

## Writes the length of genes in codon from the all_codon.snps.nonsyn file.

use strict;
use warnings;

my $usage = 'perl codon_length.pl input';

die $usage unless @ARGV;

open (IN, $ARGV[0]);
open OUT, ">$ARGV[0].len";

my @gene= ();
my $ECU = undef;

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^(ECU\d+\S+)\t(\d+)\t(\d+)/) {
		my $locus = $1;
		my $len = $3;
	
		if ($ECU =~ undef) {
		$ECU = $locus;
		print OUT "$locus\t$len\n";
		}
		elsif ($ECU eq $locus) {
		next;
		}
		elsif ($ECU ne $locus) {
		print OUT "$locus\t$len\n";
		$ECU = $locus;
		}
	}
}
close IN;
close OUT;
exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)