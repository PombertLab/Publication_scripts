#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift @ARGV){
	open IN, "<$file";
	open OUT, ">$file.parsed";
	my$ID = 'ABC';
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(CHLREDRAFT_\d{1,})/){
			my $gene = $1;
			if ($gene ne $ID){
				print OUT "$line\n";
				$ID = $gene;
			}
		}
	}
}
			