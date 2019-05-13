#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl script KEGG_file';

open IN1, "<$ARGV[0]";
open OUT, ">KEGG.hash";

my @array = ();

while (my $line = <IN1>){
	chomp $line;
	if ($line =~ /^mis:MICPUN_(\d+).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		my $stuff = "$ko\t$protein";
		push (@array, $stuff);
	}
	elsif ($line =~ /^mis:MicpuN_(\w{3}\d+).*(K\d{5})/){
		my $protein = $1;
		my $ko = $2;
		my $stuff = "$ko\t$protein";
		push (@array, $stuff);
	}
}

my @sorted = sort@array;
my $size = scalar(@sorted);
for my $count (0..$size-1){
	print OUT "$sorted[$count]\n";
}