#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.faa$//; $file =~ s/.fa$//; $file =~ s/.fasta$//;
	open OUT, ">$file.product_loci";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){print OUT "$1\n"}
	}
}