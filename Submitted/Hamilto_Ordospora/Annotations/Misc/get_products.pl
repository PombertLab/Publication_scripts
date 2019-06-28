#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'get_products.pl';
my $version = '0.1'; 

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Extracts locus_tags and their products from EMBL files
USAGE		get_products.pl -embl *.embl -out protein_list.txt
NOTE		Check the output afterwards for incomplete product names spanning two or more lines in the EMBL files

OPTIONS:
-embl		EMBL files to parse
-out		Desired output name [default: protein_list.txt]
OPTIONS
die "$usage\n" unless @ARGV;

my $output = 'protein_list.txt';
my @embl;
GetOptions(
	'out=s' => \$output,
	'embl=s@{1,}' => \@embl,
);

open OUT, ">$output";
my $key = undef;
while (my $file = shift@embl){
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
