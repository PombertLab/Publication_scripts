#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = get_products.pl *.embl";
die "\n$usage\n\n" unless @ARGV;

open OUT, ">protein_list.txt";

my $key = undef;

while (my $file = shift @ARGV){
	open IN, "<$file"; 
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^FT\s+(gene|CDS)/){$key = $1;}
		elsif ($line =~ /locus_tag=\"(.*)\"/){
			my $locus = $1;
			if ($key eq 'gene'){print OUT "$locus\thypothetical protein\n";}
		}
	}
}
