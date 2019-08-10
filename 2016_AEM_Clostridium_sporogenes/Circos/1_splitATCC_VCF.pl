#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT1, ">$file.contig1";
	open OUT2, ">$file.contig2";
	open OUT3, ">$file.contig3";
	open OUT4, ">$file.contig4";
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /ABKW02000001.1/){print OUT1 "$line\n";}
		elsif ($line =~ /ABKW02000002.1/){print OUT2 "$line\n";}
		elsif ($line =~ /ABKW02000003.1/){print OUT3 "$line\n";}
		elsif ($line =~ /ABKW02000004.1/){print OUT4 "$line\n";}
	}
}