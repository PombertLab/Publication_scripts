#!/usr/bin/perl
## Pombert Lab, IIT, 2019
## Requires tab-separated accession number / product list 
my $version = '0.2';
my $name = 'parse_UniProt_BLASTs.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a tab-delimited list of products found/not found with BLAST searches
USAGE		parse_UniProt_BLASTs.pl -b blast_output -q query.list -u uniprot_list -o parsed.tsv

OPTIONS:
-b (--blast)	BLAST tabular output (outfmt 6)
-q (--query)	List of proteins queried against UniProt
-u (--uniprot)	Tab-delimited list of UniProt accesssion numbers/products ## Created with hash_products.pl
-o (--output)	Desired output name
OPTIONS
die "$options\n" unless @ARGV;

my $blast;
my $query;
my $uniprot;
my $output;
GetOptions(
	'b|blast=s' => \$blast,
	'q|query=s' => \$query,
	'u|uniprot=s' => \$uniprot,
	'o|output=s' => \$output
);

open QUERIES, "<$query";
open PRODUCTS, "<$uniprot";
open BLAST, "$blast";
open OUT, ">$output";

my %products = (); ## Empty product hash; WARNING will take a decent chunk of RAM for 
my %hits = (); ## Keeping track of BLAST hits
my %evalues = (); 

## Filling the product hash
while (my $line = <PRODUCTS>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$products{$1}=$2;}
}

## Parsing BLAST
while (my $line = <BLAST>){
	chomp $line;
	my @columns = split("\t", $line);
	my $query = $columns[0];
	my $hit = $columns[1];
	my $evalue =$columns[10];
	if (exists $hits{$query}){next;}
	elsif ($products{$hit} =~ /uncharacterized/i){next;} ## Discarding uninformative BLAST hists
	elsif ($products{$hit} =~ /hypothetical/i){next;} ## Discarding uninformative BLAST hists
	elsif ($products{$hit} =~ /predicted protein/i){next;} ## Discarding uninformative BLAST hists
	else{
		$hits{$query}=$products{$hit};
		$evalues{$query}=$evalue;
	}
}
## Writing parsed list
print OUT "\###Protein\tE-value\tProduct\n";
while (my $line = <QUERIES>){
	chomp $line;
	if ($line =~ /^(\S+)/){
		my $query = $1;
		if (exists $hits{$query}){print OUT "$query\t$evalues{$query}\t$hits{$query}\n";}
		else{print OUT "$query\tN\/A\thypothetical protein\n";}
	}
}
