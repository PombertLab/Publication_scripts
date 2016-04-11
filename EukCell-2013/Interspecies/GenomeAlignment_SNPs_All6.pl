#!/usr/bin/perl

## Loads each sequence and put them into their own arrays with one base per line.
## This is possible because the genomes have been aligned. A single frame shift will kill this script.
## Each aligned position is compared between the three genomes, and written to the appropriate output files.
## Note: The sequence must be on one line without a fasta header. Use the fasta_to_string.pl script to convert to strings.


use strict;
use warnings;

my $usage = 'perl GenomeAlignment_SNPs_All6.pl ECI.string ECII.string ECIII.string INT.string HEL.string ROM.string';

die $usage unless @ARGV;

open (IN1, $ARGV[0]);
open (IN2, $ARGV[1]);
open (IN3, $ARGV[2]);
open (IN4, $ARGV[3]);
open (IN5, $ARGV[4]);
open (IN6, $ARGV[5]);

open INVAR, ">aligned_genomes.invar";
open SNPs, ">aligned_genomes.snps";
open GAPs, ">aligned_genomes.gaps";

my $strain1 = <IN1>;
my $strain2 = <IN2>;
my $strain3 = <IN3>;
my $strain4 = <IN4>;
my $strain5 = <IN5>;
my $strain6 = <IN6>;

my @ECI = split ('', $strain1);
my @ECII = split ('', $strain2);
my @ECIII = split ('', $strain3);
my @INT = split ('', $strain4);
my @HEL = split ('', $strain5);
my @ROM = split ('', $strain6);

my $length = scalar @ECI;
my $end = $length-1;

my $position = 1;

	foreach (0..$end) {
	
		if (($ECI[$_] eq '-') || ($ECII[$_] eq '-') || ($ECIII[$_] eq '-')|| ($INT[$_] eq '-')|| ($HEL[$_] eq '-')|| ($ROM[$_] eq '-')) {
		
		print GAPs "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\t$INT[$_]\t$HEL[$_]\t$ROM[$_]\n";
		
		$position++; 
		
		}
				
		elsif (($ECI[$_] eq $ECII[$_]) && ($ECI[$_] eq $ECIII[$_])&& ($ECI[$_] eq $INT[$_])&& ($ECI[$_] eq $HEL[$_])&& ($ECI[$_] eq $ROM[$_])) {
		
		print INVAR "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\t$INT[$_]\t$HEL[$_]\t$ROM[$_]\n";
		
		$position++; 
		
		}
		
		else {
		
		print SNPs "$position\t$ECI[$_]\t$ECII[$_]\t$ECIII[$_]\t$INT[$_]\t$HEL[$_]\t$ROM[$_]\n";
		
		$position++; 
		
		}
	}

close IN1;
close IN2;
close IN3;
close IN4;
close IN5;
close IN6;
close INVAR;
close SNPs;
close GAPs;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)