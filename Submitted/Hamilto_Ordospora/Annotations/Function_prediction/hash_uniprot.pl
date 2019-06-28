#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'hash_uniprot.pl';
my $version = '0.1';

use strict; use warnings;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Creates tab-delimited lists of products from the trEMBL and SwissProt FASTA files
USAGE		hash_uniprot.pl uniprot_sprot.fasta uniprot_trembl.fasta
OPTIONS
die "$usage\n" unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.fasta$//;
	open OUT, ">$file.list";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){
			print OUT "$1\t$2\n";
		}
	}
}
