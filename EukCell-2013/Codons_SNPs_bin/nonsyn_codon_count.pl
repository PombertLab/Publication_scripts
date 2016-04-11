#!/usr/bin/perl

## Use this script on all_codon.snps.nonsyn to count the number of non-synonymous changes per gene.
## Works also with the all_codon.snps.syn as input to count synonymous changes.
## The output is out of order from the hash, just use Excel of LibreOffice to sort it out.

use strict;
use warnings;

my $usage = 'perl nonsyn_codon_count.pl input';

die $usage unless @ARGV;

open (IN, $ARGV[0]);
open OUT, ">$ARGV[0].count";

my @gene= ();
my %seen = ();

while (my $line = <IN>) {
	chomp $line;
	if ($line =~ /^(ECU\d+\S+)\t(\d+)/) {
	push (@gene, $1);
	}
}

foreach my $ECU (@gene){
  $seen{$ECU}++
}

foreach my $ECU (keys %seen){
  print OUT "$ECU\t$seen{$ECU}\n"
}

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)