#!/usr/bin/perl

use strict;
use warnings;

open IN, "<$ARGV[0]";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^(\S+)/){
		my $contig = $1;
		open OUT, ">>$contig.gff3";
		print OUT "$line\n";
		close OUT;
	}
}