#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'split_CSV.pl';
my $version = 0.1;

use strict; use warnings;

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Splits single copy and multicopy genes from Orthofinder into separate CSV files
USAGE		split_CSV.pl Orthogroups.csv
OPTIONS
die "$options\n" unless @ARGV;

system "dos2unix $ARGV[0]"; ## Removing DOS format SNAFUs...
open IN, "<$ARGV[0]";
open ORT, ">single_copy.csv";
open PAR, ">multi_copy.csv";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^\t/){print ORT "$line\n"; print PAR "$line\n";}
	elsif ($line =~ /,/){print PAR "$line\n";}
	else {print ORT "$line\n";}
}
