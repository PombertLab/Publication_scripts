#!/usr/bin/perl
## Pombert Lab, 2017
my $version = '0.1';
my $name = 'hash_products.pl';

use strict; use warnings;

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates lists of products from trEMBL and SwissProt fasta headers
USAGE		hash_products.pl *.fasta
OPTIONS
die "$options\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fasta$//;
	open OUT, ">$file.products.hash";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){
			print OUT "$1\t$2\n";
		}
	}
}
