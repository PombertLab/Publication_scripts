#!/usr/bin/perl
## Pombert lab 2018
my $version = 0.1;
my $name = 'genes_expressed.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Uses the output of RNAseq mapping, samtools, and EMBL annotations to identify which genes are expressed
USAGE		genes_expressed.pl -c PASS.coverage -p Verified_products_ALL.list -e *.embl -o expressed_proteins.tsv

-c	Coverage file obtained with Samtools # samtools depth -aa PASS_sorted.bam > PASS.coverage
-p	Product list (tab-delimited)
-e	Annotations (in EMBL format)
-o	Output file in TSV format [Default = expressed_proteins.tsv]
OPTIONS
die "$usage\n" unless @ARGV;

my $coverage;
my $products;
my @embl;
my $output = 'expressed_proteins.tsv';
GetOptions(
	'c=s'=>\$coverage,
	'p=s'=>\$products,
	'o=s'=>\$output,
	'e=s@{1,}'=>\@embl
);


## Creating coverage database
open COV, "<$coverage" or die "Can't open $coverage!\n";
my %coverage;
system "echo Working on coverage database...";
while (my $line = <COV>){
	chomp $line;
	if ($line =~ /^(\w+)\t(\d+)\t(\d+)/){
		my $contig = $1; #print "$contig\t";
		my $position = $2; #print "$position\t";
		my $cov = $3; #print "$cov\n";
		$coverage{$contig}{$position} = $cov;
	}
}
## Creating products database
open PROD, "<$products" or die "Can't open $products!\n";
my %products;
while (my $line = <PROD>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif($line =~ /^(\w+)\t(.*)$/){$products{$1}=$2;}
}

## Working through CDS annotations; disregarding tRNA and rRNAs...
system "echo Working on EMBL files...";
open OUT, ">$output";
print OUT "Protein\tProduct\tSum\tLength (in nt)\tAverage RNAseq depth\n";
my $sum; my $size; my $average; my $cg; my $locus_tag;
while (my $embl = shift@embl){
	open EMBL, "<$embl";
	if ($embl =~ /^(\w+)/){$cg = $1;}
	while (my $line = <EMBL>){
		chomp $line;
		$sum = 0; $size = 0; $average = 0;
		if ($line =~ /^FT\s+\/locus_tag="(.*)"/){$locus_tag = $1;}
		elsif ($line =~ /^FT \s+CDS\s+.*\((\S+)\)/){depth();}
		elsif ($line =~ /^FT \s+CDS\s+(.*)$/){depth();}
	}
}

## Subroutines
sub depth {
	my @nucleotides = split(",", $1);
		while (my $nt = shift@nucleotides){
			if ($nt =~ /^(\d+)\.\.(\d+)/){
				my $start = $1; my $end = $2;
				for my $pos ($start..$end){
					$sum += "$coverage{$cg}{$pos}";
					$size++;
				}
			} 
		}
	$average = sprintf ("%.2f", $sum/$size); # unless $size == 0;
	print OUT "$locus_tag\t$products{$locus_tag}\t$sum\t$size\t$average\n";
	$sum = 0; $size = 0; $average = 0;
}