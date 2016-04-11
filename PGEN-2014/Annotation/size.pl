#!/usr/bin/perl

use strict;
use warnings;

open OUT, ">Contig_sizes.txt";

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.string//;
	while (my $line = <IN>){
		chomp $line;
		my $size=length($line);
		print OUT "$file\t$size\n";
	}
}
