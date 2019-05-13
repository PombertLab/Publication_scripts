#!/usr/bin/perl

## Pombert Lab, IIT 2014
## Parses InterProScan for GOs and other relevant information

use strict;
use warnings;

my $usage = 'perl filter_InterProScan.pl *.interpro.tsv'; die $usage unless@ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.tsv$//;
	open GO, ">$file.GOs";
	
	my %hash = ();

	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /GO/){
			my @array = split("\t", $line);
			my $locus = $array[0]; # $ARGV[0]; Protein Accession (e.g. P51587)
			# Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
			# Sequence Length (e.g. 3418)
			# Analysis (e.g. Pfam / PRINTS / Gene3D)
			# Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
			# Signature Description (e.g. BRCA2 repeat profile)
			# Start location
			# Stop location
			# Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
			# Status - is the status of the match (T: true)
			# Date - is the date of the run
			# $ARGV[11] (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
			# $ARGV[12] (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
			# $ARGV[13]; (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
			# $ARGV[14]; (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
			$array[13] =~ s/\|/ /g;
			push (@{$hash{$locus}}, "$array[13]");
		}
	}
	
	for my $key ( sort keys %hash ) {
		print GO "$key\t@{$hash{$key}}\n";
	}
}
			


