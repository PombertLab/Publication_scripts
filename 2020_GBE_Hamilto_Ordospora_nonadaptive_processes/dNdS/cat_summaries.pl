#!/usr/bin/perl
## Pombert Lab, 2018
my $name = 'cat_summaries.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Quick script to generate a summary table from the SNAP summary files
USAGE		cat_summaries.pl -p products.txt -spp species.txt -sum summary.* -o summaries_table.tsv

OPTIONS:
-p	## Tab-delimited list of locus tags and their products [defailt: products.txt]
-spp	## Tab-delimited list of locus tag prefixes and corresponding species [default: species.txt]
-sum	## Summary output files created by run_SNAP.pl/SNAP_mod.pl
-o	## Desired output name for the tab-delimited table
OPTIONS
die "$usage\n" unless @ARGV;

my $products = 'products.txt';
my $species = 'species.txt';
my @summaries;
my $output;
GetOptions(
	'p=s' => \$products,
	'spp=s' => \$species,
	'sum=s@{1,}' => \@summaries,
	'o=s' => \$output
);

## Products database
open PD, "<$products";
my %prod;
while (my $line = <PD>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$prod{$1} = $2;}
}
## Species database
open SPP, "<$species" or die "Can't open file $species\n";
my %species;
while (my $line = <SPP>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$species{$1} = $2;}
}

## Working on summaries
open OUT, ">$output";
print OUT "Orthogroup\tSpecies_1\tseq1\tProduct(seq1)\tSpecies_2\tseq2\tProduct(seq2)\tds\tdn\tds\/dn\tdn\/ds\n";
while (my $summary = shift@summaries){
	open SUM, "<$summary";
	my $og = $summary; $og =~ s/^summary.//; $og =~ s/.NT$//; $og =~ s/.AA$//;
	while (my $line = <SUM>){
		chomp $line;
		if ($line =~ /^Compare/){next;}
		elsif ($line =~ /^Averages/){next;}
		else{
			my @col = split(/\s+/, $line);
			# $col[2] = sequence 1
			# $col[3] = sequence 2
			# $col[10] = ds; $col[11] = dn; $col[12] = ds/dn;
			my $sp1; my $sp2;
			if ($col[2] =~ /^(\S+_)/){$sp1 = $1;}
			if ($col[3] =~ /^(\S+_)/){$sp2 = $1;}
			my $dsdn;
			if (($col[10] eq 'NA') || ($col[11] eq 'NA')){$dsdn = 'NA';}
			elsif (($col[10] <= 0) || ($col[11] <= 0)){$dsdn = 'NA';}
			elsif (($col[10] > 0) && ($col[11] > 0)){$dsdn = ($col[11]/$col[10]); $dsdn = sprintf("%.4f", $dsdn);}
			print OUT "$og\t$species{$sp1}\t$col[2]\t$prod{$col[2]}\t$species{$sp2}\t$col[3]\t$prod{$col[3]}\t$col[10]\t$col[11]\t$col[12]\t$dsdn\n";
		}
	}
}