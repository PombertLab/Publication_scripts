#!/usr/bin/perl

## This script generates a bin file to be used with SNPs_slide.pl
## Loads each sequence ECI, ECII and ECIII and put them into their own arrays with one base per line.
## This is possible because the genomes have been aligned. A single frame shift will kill this script.
## Each aligned position is compared between the three genomes, and written 0 or 1 to the bin file. 0 = invariable, 1 = SNP.
## The string from the bin file is to be used with the script for the sliding windows.

use strict;
use warnings;

my $usage = 'perl Binary_SNPs.pl ECI.string ECII.string ECIII.string';

die $usage unless @ARGV;

open (IN1, $ARGV[0]);
open (IN2, $ARGV[1]);
open (IN3, $ARGV[2]);

open BIN, ">aligned_genomes.bin"; ## Replace with desired output name

my $strain1 = <IN1>;
my $strain2 = <IN2>;
my $strain3 = <IN3>;

my @ECI = split ('', $strain1);
my @ECII = split ('', $strain2);
my @ECIII = split ('', $strain3);

my $length = scalar @ECI;
my $end = $length-1;

	foreach (0..$end) {
				
		if (($ECI[$_] eq $ECII[$_]) && ($ECI[$_] eq $ECIII[$_])) {
			print BIN "0";
		}
		
		else {
			print BIN "1";
		}
	}

close IN1;
close IN2;
close IN3;
close BIN;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)
