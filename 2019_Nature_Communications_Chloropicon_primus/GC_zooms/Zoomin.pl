#!/usr/bin/perl
## Pombert Lab, 2019
my $name = 'Zoomin.pl';
my $version = 0.2;

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Calculates the GC content of a locus and returns the genes contained within the locus
USAGE		Zoomin.pl -f chromosome_XX.fsa -s start -e end -l Verified_products_ALL.list -embl chromosome_XX.embl 

EXAMPLE		Zoomin.pl -f Chromosome_03.fsa -s 1337500 -e 1369999 -l Verified_products_ALL.list -embl Chromosome_03.embl 

OPTIONS:
-f (--fasta)	FASTA file of the chromosome or contig
-s (--start)	Start of the region to look at
-e (--end)	End of the region to look at
-l (--list)	List of locus_tags and their products
-embl		Annotations in EMBL format
OPTIONS
die "$options\n" unless @ARGV;

my $products;
my $start;
my $end;
my $fasta;
my $embl;

GetOptions(
	'l|list=s' => \$products,
	's|start=i' => \$start,
	'e|end=i' => \$end,
	'f|fasta=s' => \$fasta,
	'embl=s' => \$embl
);

## Creating database of products
my %products;
open PROD, "<$products" or die "Can\'t open products file $products...\n";
while (my $line = <PROD>){
	chomp $line;
	if ($line =~ /^#/){next;}
	else{
		my @array = split("\t", $line);
		$products{$array[0]} = $array[1];
	}
}

## Calculating GC content in the zoomed in region
open FASTA, "<$fasta" or die "Can\'t open fasta file $fasta...\n";
my $sequence;
while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>/){next;}
	else {$sequence .= $line;}
}
my $sub = substr($sequence, $start - 1, $end - $start - 1);
my $size = length$sub;
my $GC = 0;
my @nt = unpack("(A1)*", $sub);
while (my $base = shift@nt){if ($base =~ /[GCgc]/){$GC++;}}
my $percent = sprintf("%.1f", ($GC/$size)*100);

## Extracting genes located in the zoomed in region (genes across the boundaries of the regions will not be printed)  
my $gene_start; my $gene_end;
my $locus; my $key;
my $strand;
open EMBL, "<$embl" or die "Can\'t open EMBL file $embl...\n";
$embl =~ s/.embl$//;
open OUT, ">$embl.zoomin";
print OUT "## Zoom in on $embl from position $start to position $end; Length $size, GC% $percent\n";
print OUT "Locus tag\tstrand\tStart\tEnd\tProduct\n";
while (my $line = <EMBL>){
	chomp $line;
	if ($line =~ /^FT\s+(gene)\s+(\d+)\.\.(\d+)/){$key = $1; $gene_start = $2; $gene_end = $3; $strand = 'plus';}
	elsif ($line =~ /^FT\s+(gene)\s+complement\((\d+)\.\.(\d+)/){$key = $1; $gene_start = $3; $gene_end = $2; $strand = 'minus';}
	elsif ($line =~ /^FT\s+(CDS|rRNA|tRNA)/){$key = $1;}
	if (($key eq 'gene') && ($line =~ /^FT\s+\/locus_tag=\"(\w+)\"/)) {
		$locus = $1;
		if ($strand eq 'plus'){if (($gene_start >= $start) && ($gene_end <= $end)) {print OUT "$locus\t$strand\t$gene_start\t$gene_end\t$products{$locus}\n";}}
		if ($strand eq 'minus'){if (($gene_end >= $start) && ($gene_start <= $end)) {print OUT "$locus\t$strand\t$gene_start\t$gene_end\t$products{$locus}\n";}}
	}
}