#!/usr/bin/perl

## 	Hash = ECU_products_list.txt (a tab delimited gene & product list, i.e. =~ /gene\tproduct\n/ )
##	Query = output of sort_SNP_metrics.pl

use strict;
use warnings;

my $usage = 'perl bin_products.pl hash query';

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
open OUT, ">$ARGV[1].out";
print OUT "Gene\tlength\tSNPs\tProduct\n";

	while (my $query = <QUERY>) {
		chomp $query;
		
		if ($query =~ /^(ECU\d+_\S+)\t(\d+)\t(\d+)\t\d.*/) {
			my $gene = $1;
			my $len = $2;
			my $snps = $3;
			print OUT "$gene\t$len\t$snps\t$hash{$gene}\n";
		}
	}


close HASH;
close QUERY;
close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)