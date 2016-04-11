#!/usr/bin/perl

use strict;
use warnings;

while (my $input = shift @ARGV){
	
	open IN, "<$input";
	open OUT, ">$input.filtered";
	my $locus = 'ABC';
	
	while (my $line = <IN>){
		
		chomp $line;
		
		if (($line =~ /^(\d{1,})\t/) && ($locus eq 'ABC')){
			my $id = $1;
			$locus = $id;
			print OUT "$line\n";
		}
		elsif (($line =~ /^(\d{1,})\t/) && ($locus ne 'ABC')){
			my $id = $1;
				if ($id eq $locus){
					next;
				}
				elsif ($id ne $locus){
					print OUT "$line\n";
					$locus = $id;
				}
		}
	}
}