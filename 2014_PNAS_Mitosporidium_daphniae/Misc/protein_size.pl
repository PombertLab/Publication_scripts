#!/usr/bin/perl

use strict;
use warnings;

open OUT, ">protein_sizes.txt";

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fsa$//;
	my $protein = undef;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(DI09_\d+p\d+)/ ){
			$protein=$1;
			print OUT "$protein\t";
		}
		else{
			my $size = length($line);
			print OUT "$size\n";
		}
	}
}