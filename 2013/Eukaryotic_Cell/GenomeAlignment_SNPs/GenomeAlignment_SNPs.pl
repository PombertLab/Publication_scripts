#!/usr/bin/perl

## Loads each sequence ECI, ECII and ECIII and put them into their own arrays with one base per line.
## This is possible because the genomes have been aligned. A single frame shift will kill this script.
## Each aligned position is compared between the three genomes, and written to the appropriate output files.
## Note: The sequence must be on one line without a fasta header. Use the fasta_to_string.pl script to convert to strings.


use strict;
use warnings;

my $usage = 'perl GenomeAlignment_SNPs.pl ECI.string ECII.string ECIII.string';

die $usage unless @ARGV;

open (IN1, $ARGV[0]);
open (IN2, $ARGV[1]);
open (IN3, $ARGV[2]);

open INVAR, ">aligned_genomes.invar";
open SNPs, ">aligned_genomes.snps";
open GAPs, ">aligned_genomes.gaps";

my $strain1 = <IN1>;
my $strain2 = <IN2>;
my $strain3 = <IN3>;

my @ECI = split ('', $strain1);
my @ECII = split ('', $strain2);
my @ECIII = split ('', $strain3);

my $length = scalar @ECI;
my $end = $length-1;

my $position = 1;

	foreach (0..$end) {
	
		if (($ECI[$_] eq '-') || ($ECII[$_] eq '-') || ($ECIII[$_] eq '-')) {
		
		print GAPs "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
		
		$position++; 
		
		}
				
		elsif (($ECI[$_] eq $ECII[$_]) && ($ECI[$_] eq $ECIII[$_])) {
		
		print INVAR "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
		
		$position++; 
		
		}
		
		else {
		
		print SNPs "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\n";
		
		$position++; 
		
		}
	}

close IN1;
close IN2;
close IN3;
close INVAR;
close SNPs;
close GAPs;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)