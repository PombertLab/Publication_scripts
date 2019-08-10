#!/usr/bin/perl

use strict;
use warnings;

## use .len as the hash (gene length in codons), the .count file as query

my $usage = 'perl nonsyn_codon_average.pl hash query';

die $usage unless @ARGV;

open (HASH, $ARGV[0]);
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
print OUT "Gene\tlength\tSNPs\tRatio\n";

	while (my $query = <QUERY>) {
		chomp $query;
		
		if ($query =~ /^(ECU\d+_\S+)\t(\d+)/) {
			my $gene = $1;
			my $snps = $2;
			my $av = ($snps/$hash{$gene});
			print OUT "$gene\t$hash{$gene}\t$snps\t$av\n";
		}
	}


close HASH;
close QUERY;
close OUT;

exit;

## Written by JF Pombert, Canadian Institute for Advanced Research
## University of British Columbia (2012)