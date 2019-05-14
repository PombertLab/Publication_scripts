#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.1';
my $name = 'split_VCF.pl';

use strict; use warnings;

# Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Splits VCF file(s) per contig/chromosome, excluding positions masked (N) by RepeatMasker in the reference.
USAGE		sort_VCF.pl *.vcf
OPTIONS
die "$options\n" unless @ARGV;

## Parsing VCF files
while (my $vcf = shift@ARGV){
	open VCF, "<$vcf";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @columns = split("\t", $line);
			unless ($columns[3] eq 'N'){ ## Checking if position has been masked by RepeatMasker, if not continue...
				open OUT, ">>$columns[0].vcf";
				print OUT "$line\n";
				close OUT;
			}
		}
	}
}