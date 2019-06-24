#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'USAGE = txid.pl input.fasta txid_number output_list';
die "\n$usage\n\n" unless @ARGV;

open IN, "<$ARGV[0]";
my $txid = $ARGV[1];
open OUT, ">$ARGV[2]";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>(\S+)/){
		print OUT "$1\t$txid\n";
	}
}
