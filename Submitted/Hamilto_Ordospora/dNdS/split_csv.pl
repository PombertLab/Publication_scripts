#!/usr/bin/perl
## splits single copy genes and mutlicopy genes into their own files

use strict;
use warnings;

die "USAGE = split_csv.pl Orthogroups.csv\n" unless @ARGV;

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
