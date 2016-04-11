#!/usr/bin/perl

## Generate a records file to extract desired proteins with faSomeRecords [available from UCSC]

use warnings;
use strict;

open(IN, $ARGV[0]);
open OUT, ">otauri.records";

while(my $line = <IN>){
	chomp $line;
	if($line =~ /^(Ot\d{1,}g\d{1,})/){
		my $gene = $1;
		print OUT "$gene\n";
	}
}