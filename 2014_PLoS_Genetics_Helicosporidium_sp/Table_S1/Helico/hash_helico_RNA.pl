#!usr/bin/perl

use strict;
use warnings;

while (my $file = shift @ARGV){
	open IN, "<$file";
	open OUT, ">$file.hash";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(CHLREDRAFT_\d+)\t(a\S+)/){
			my $chlamy = $1;
			my $helico = $2;
			print OUT "$chlamy\t$helico\n";
		}
	}
}