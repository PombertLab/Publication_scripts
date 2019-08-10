#!/usr/bin/perl

## Requires BLAST tabular output (outfmt 6)
## Requires tab-separated accession number / product list 

use strict;
use warnings;

my $usage = 'USAGE = perl parse_UniProt_BLASTs.pl query_list Uniprot_product_list blast_output';
die "\n$usage\n
query_list= list of proteins queried against UniProt
Uniprot_product_list= tab-delimited UniProt accesssion number/product list
blast_output = BLAST tabular output (outfmt 6)\n
" unless (scalar@ARGV==3);

open QUERIES, "<$ARGV[0]";
open PRODUCTS, "<$ARGV[1]";
open BLAST, "<$ARGV[2]";
open OUT, ">$ARGV[2].parsed";

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
	if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/){
		my $query = $1;
		my $hit = $2;
		my $evalue =$11;
		if (exists $hits{$query}){next;}
		#elsif ($products{$hit} =~ /uncharacterized/i){next;} ## Discarding uninformative BLAST hists
		#elsif ($products{$hit} =~ /hypothetical/i){next;} ## Discarding uninformative BLAST hists
		#elsif ($products{$hit} =~ /predicted protein/i){next;} ## Discarding uninformative BLAST hists
		else{
			$hits{$query}=$products{$hit};
			$evalues{$query}=$evalue;
		}
	}
}
## Writing parsed list
print OUT "\###Protein\tE-value\tProduct\n";
while (my $line = <QUERIES>){
	chomp $line;
	if ($line =~ /^(\S+)/){
		my $query = $1;
		if (exists $hits{$query}){
			print OUT "$query\t$evalues{$query}\t$hits{$query}\n";
		}
		else{
			print OUT "$query\tN\/A\thypothetical protein\n";
		}
	}
}
