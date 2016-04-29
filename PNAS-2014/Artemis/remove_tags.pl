#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Removes locus tags from EMBL files. Useful to relabel them.

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.2";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /FT                   \/locus_tag=/){next;}
		else {print OUT "$line\n";}
	}
}