#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'get_protein_loci.pl';
my $version = 0.1;

use strict; use warnings;

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generate lists of loci (one per line) from FASTA headers
EXAMPLE		get_protein_loci.pl *.fasta
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