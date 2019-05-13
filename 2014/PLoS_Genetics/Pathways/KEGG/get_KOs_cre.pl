#!/usr/bin/perl

use strict;
use warnings;

open (IN, $ARGV[0]);
open OUT, ">$ARGV[0].KOs";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /(CHLREDRAFT_\d{1,})/){
		my $gene = $1;
		print OUT "$gene\t";
	}
	elsif ($line =~/(\[KO:.\S+)/){
		my $EC = $1;
		print OUT "$EC\n";
	}
}