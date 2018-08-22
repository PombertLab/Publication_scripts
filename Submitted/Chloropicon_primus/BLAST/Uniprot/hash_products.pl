#!/usr/bin/perl
## Get products list from trEMBL and SwissProt fasta headers

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fasta$//;
	open OUT, ">$file.products.hash";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){
			print OUT "$1\t$2\n";
		}
	}
}
