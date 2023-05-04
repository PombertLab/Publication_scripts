#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2022
my $name = 'get_prod_from_gbk.pl';
my $version = '0.1';
my $updated = '2022-10-21';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage =<<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates a tab-delimited locus_tag => product from the genes annotated
		in a GenBank file.
		
COMMAND		${name} -g file.gbk -t products.tsv

OPTIONS:
-g (--gbk)	GenBank file (.gb|gbk|gbff) to parse
-t (--tsv)	Desired output TSV file 
USAGE
die "\n$usage\n" unless @ARGV;

my $gbk;
my $tsv;
GetOptions(
	'g|gbk=s' => \$gbk,
	't|tsv=s' => \$tsv
);

open GBK, "<", $gbk or die "Can't open $gbk: $!\n";
open TSV, ">", $tsv or die "Can't open $tsv: $!\n";

## Creating a single string from the GBK file
my $string;
while (my $line = <GBK>){
	chomp $line;
	$string .= $line;
}

## Spliting string per gene entry
my @data = split(/     gene            /, $string);

## Print out locus_tag + product data for each gene entry
while (my $chunk = shift@data){
	if ($chunk =~ /locus_tag=\"(\w+)\".*\/product=\"(.*?)\"/){
		my $locus = $1;
		my $product = $2;
		$product =~ s/\s{2,}/ /g;
		print TSV "$locus\t$product\n";
	} 
}