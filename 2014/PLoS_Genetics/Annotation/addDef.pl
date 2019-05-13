#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.added";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>/){
			print OUT "$line".'[host=Simulium jonesi][country=USA:Florida][collection_date=1998][isolation_source=black fly larvae]'."\n";
		}else{
			print OUT "$line\n";
		}
	}
}