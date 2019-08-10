#!/usr/bin/perl

## Generate a records file to extract desired proteins with faSomeRecords [available from UCSC]

use warnings;
use strict;

open(IN, $ARGV[0]);
open OUT, ">volvox.records";

while(my $line = <IN>){
	chomp $line;
	if($line =~ /^(VOLCADRAFT_\d{1,})/){
		my $gene = $1;
		print OUT "$gene\n";
	}
}