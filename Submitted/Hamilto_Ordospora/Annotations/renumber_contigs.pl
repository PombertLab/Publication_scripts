#!/usr/bin/perl

use strict;
use warnings;

my $usage = "rename.pl fasta_file output_file";
die "\n$usage\n\n" unless @ARGV;

open IN, "<$ARGV[0]";
open OUT, ">$ARGV[1]";

my $count = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>/){
	 	$count++;
		my $num = sprintf("%04d", $count);
		print OUT ">contig-$num\n";
	}
	else {
		print OUT "$line\n";
	}
}
