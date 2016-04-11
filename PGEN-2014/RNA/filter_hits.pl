#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift @ARGV){
	open IN, "<$file";
	open OUT, ">$file.filtered";
	my $ID = 'undef';
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(a\d+_\d+)/){
			my $RNA = $1;
			if ($RNA eq $ID){
				next;
			}
			elsif ($RNA ne $ID){
				print OUT "$line\n";
				$ID = $RNA;
			}
		}
	}
}