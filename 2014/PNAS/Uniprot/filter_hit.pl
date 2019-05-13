#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.filtered";
	my $protein = undef;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(DI09_\d+p\d+)\ttr\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\S+\t\d+e-0[1234]/){next;}
		elsif ($line =~ /^(DI09_\d+p\d+)/){
			my $query = $1;
			if (!defined$protein){print OUT "$line\n";$protein=$query;}
			elsif ($query eq $protein){next;}
			else {print OUT "$line\n"; $protein=$query;}
		}
	}
}