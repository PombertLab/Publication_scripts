#!/usr/bin/perl

use strict;
use warnings;

die "\nUSAGE = low_gc.pl maxGC genome.gc\n\n" unless @ARGV;

my $gc = $ARGV[0];
my $file = $ARGV[1];

open IN, "<$file";
open OUT, ">$file.$gc.low_gc";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	else{
		my @cols = split(" ", $line);
		if ($cols[3] <= $gc){print OUT "$line\n";}
	}
}
