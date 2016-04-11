#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift @ARGV){
	open IN, "<$file";
	open OUT, ">$file.oneliner";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>/){
			print OUT "\n$line\n";
		}
		else {
			print OUT "$line";
		}
	}
}