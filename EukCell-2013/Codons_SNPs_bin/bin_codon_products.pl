#!/usr/bin/perl

## 	Hash = ECU_products_list.txt (a tab delimited gene & product list, i.e. =~ /gene\tproduct\n/ )
##	Query = output of nonsyn_codon_average.pl

use strict;
use warnings;

my $usage = 'perl bin_codon_products.pl hash query';

die $usage unless @ARGV;

open (HASH, $ARGV[0]); ## The hash input is a gene & product tab-delimited list, one per line,.
my %hash = ();

#### Fill the hash with the genes (keys) and products (values).

	while (my $line = <HASH>) {
		chomp $line;
		
		if ($line =~ /^(ECU\d+_\S+)\t(.*)/) {
			my $key = $1;
			my $value = $2;
			$hash{"$key"} = "$value";
		}
	}

#### Work on the query file.

open (QUERY, $ARGV[1]);
open LOW, ">$ARGV[1].low";
open MID, ">$ARGV[1].mid";
open HIGH, ">$ARGV[1].high";

print LOW "Gene\tlength\tSNPs\tRatio\tProduct\n";
print MID "Gene\tlength\tSNPs\tRatio\tProduct\n";
print HIGH "Gene\tlength\tSNPs\tRatio\tProduct\n";

my $high = 0.01;	## Change values to the desired ratios
my $mid = 0.001;	## Change values to the desired ratios
my $low = 0;		## Change values to the desired ratios

	while (my $query = <QUERY>) {
		chomp $query;
		
		if ($query =~ /^(ECU\d+_\S+)\t(\d+)\t(\d+)\t(\d.\d+)/) {
			my $gene = $1;
			my $len = $2;
			my $snps = $3;
			my $ratio = $4;
			
			if ($ratio >= $high) {
			print HIGH "$gene\t$len\t$snps\t$ratio\t$hash{$gene}\n";
			}
			elsif (($ratio >= $mid) && ($ratio < $high)) {
			print MID "$gene\t$len\t$snps\t$ratio\t$hash{$gene}\n";
			}
			elsif (($ratio >= $low) && ($ratio < $mid)) {
			print LOW "$gene\t$len\t$snps\t$ratio\t$hash{$gene}\n";
			}
		}
	}


close HASH;
close QUERY;
close LOW;
close MID;
close HIGH;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)