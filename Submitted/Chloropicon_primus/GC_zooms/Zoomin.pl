#!/usr/bin/perl
## Pombert Lab, 2017
## v0.1

use strict;
use warnings;

die "\nUSAGE = Zoomin.pl Verified_products_ALL.list start end chromosome_XX.embl chromosome_XX.fsa\n\n" unless @ARGV;

my $products = $ARGV[0];
my $start = $ARGV[1];
my $end = $ARGV[2];
my $chromo = $ARGV[3];
my $fasta = $ARGV[4];

## Creating database of products
my %products;
open PROD, "<$products";
while (my $line = <PROD>){
	chomp $line;
	if ($line =~ /^#/){next;}
	else{
		my @array = split("\t", $line);
		$products{$array[0]} = $array[1];
	}
}

## Calculating GC content in the zoomed in region
open FASTA, "<$fasta";
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
open EMBL, "<$chromo";
$chromo =~ s/.embl$//;
open OUT, ">$chromo.zoomin";
print OUT "## Zoom in on $chromo from position $start to position $end; Length $size, GC% $percent\n";
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