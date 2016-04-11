#!/usr/bin/perl

## Generate a records file to extract desired proteins with faSomeRecords [available from UCSC]

use warnings;
use strict;

open(IN, $ARGV[0]);
open OUT, ">lucima.records";

while(my $line = <IN>){
	chomp $line;
	if($line =~ /^(OSTLU_\d{1,})/){
		my $gene = $1;
		print OUT "$gene\n";
	}
}