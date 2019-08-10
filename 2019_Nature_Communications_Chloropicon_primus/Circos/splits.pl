#!/usr/bin/perl
## Pombert Lab, 2017
my $name = 'splits.pl';
my $version = 0.1;

use strict; use warnings;

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Splits genome.list with KEGG orthologs in 3 files; plus strand, minus strand, and no KEGG ortholog 
USAGE		splits.pl genome.list
OPTIONS
die "$options\n" unless @ARGV;

open IN, "<$ARGV[0]";
open PLUS, ">genome.plus";
open MINUS, ">genome.minus";
open NA, ">genome.no_ko";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^\w+\t(\w+)\t(\d+)\t(\d+)\t(minus|plus)\t(K\d+|NA)/){
		my $locus = $1;
		my $start = $2;
		my $end = $3;
		my $strand = $4;
		my $ko = $5;
		if ($strand eq 'plus') {print PLUS "$locus $start $end fill_color=blue\n";}
		if ($strand eq 'minus') {print MINUS "$locus $start $end fill_color=red\n";}
		if ($ko eq 'NA'){print NA "$locus $start $end fill_color=grey\n";}
	}
}