#!/usr/bin/perl
## Pombert Lab, 2017
my $version = '0.1';
my $name = 'get_protein_loci.pl';

use strict; use warnings;

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Extract locus tags (protein loci) from fasta headers
USAGE		get_protein_loci.pl *.fasta
OPTIONS
die "$options\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.faa$//; $file =~ s/.fa$//; $file =~ s/.fasta$//;
	open OUT, ">$file.product_loci";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){print OUT "$1\n"}
	}
}