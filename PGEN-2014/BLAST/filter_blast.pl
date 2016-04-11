#!/usr/bin/perl

use strict;
use warnings;

while (my $input = shift @ARGV){
	
	open IN, "<$input";
	open OUT, ">$input.filtered";
	my $locus = 'ABC';
	
	while (my $line = <IN>){
		
		chomp $line;
		
		if (($line =~ /^(H632_c\d{1,}\.\d{1,})/) && ($locus eq 'ABC')){
			my $id = $1;
			$locus = $id;
			print OUT "$id\n";
		}
		elsif (($line =~ /^(H632_c\d{1,}\.\d{1,})/) && ($locus ne 'ABC')){
			my $id = $1;
				if ($id eq $locus){
					next;
				}
				elsif ($id ne $locus){
					print OUT "$id\n";
					$locus = $id;
				}
		}
	}
}