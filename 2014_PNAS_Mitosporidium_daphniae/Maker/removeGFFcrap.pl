#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.trimmed";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^contig-\d+\t(blastx|snap_masked|augustus_masked|protein2genome)/){
			next;
		}else{print OUT "$line\n";}
	}
}