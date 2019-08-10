#!/usr/bin/perl

## 	Hash = ECUs_to_UniProt.tab (a tab delimited ECU & UniProt list, i.e. =~ /gene\tUniProt\n/ )

use strict;
use warnings;

my $usage = 'perl get_UniProt_descriptions.pl ECUs_to_UniProt.tab query';

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
open OUT, ">$ARGV[1].uni";
my $noUni = 0;
			
	while (my $query = <QUERY>) {
		chomp $query;
		
		if ($query =~ /^(ECU\d+_\S+)\t/) {
			my $gene = $1;

			
			if ($hash{$gene} eq defined) {
				$noUni++;
				next;
			}
			else {
			print OUT "$hash{$gene}\n";
			}
		}
	}

print OUT "\n$noUni genes without UniProt descriptors\n";


close HASH;
close QUERY;
close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)