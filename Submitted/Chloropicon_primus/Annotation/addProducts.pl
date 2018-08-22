#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'USAGE = perl addProducts.pl product_list fasta_file';
die "$usage\n" unless (scalar(@ARGV == 2));

open PRODUCTS, "<$ARGV[0]";
open FASTA, "<$ARGV[1]";
open OUT, ">$ARGV[1].products";

my %products = ();

while (my $line = <PRODUCTS>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^(\S+)\t(.*)$/){$products{$1} = $2;}
}

while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){
		my $locus = $1;
		print OUT ">$locus [product=$products{$locus}]\n";
	}
	else {print OUT "$line\n";}
}