#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.added";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>/){
			print OUT "$line".'[uncultured]'."\n";
		}else{
			print OUT "$line\n";
		}
	}
}