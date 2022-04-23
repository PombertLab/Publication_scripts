#!/usr/bin/perl
## Pombert Lab 2019
## Uses the output of RNAseq mapping, samtools, and GFF annotations to identify which genes are expressed
my $version = "0.5";
my $name = 'genes_expressed.pl';
my $updated = "2021-09-17";

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Calculates average RNAseq depth per genes from samtools depth files and NCBI gff annotations files

USAGE		${name} \\
		  -c PASS.coverage \\
		  -o expressed_proteins.tsv \\
		  -g *.gff

OPTIONS:
-c	Coverage file obtained with Samtools # samtools depth -aa PASS_sorted.bam > PASS.coverage
-o	Output file in TSV format [Default = expressed_genes.tsv]
-g	Annotations (in GFF format)

OPTIONS

die "$usage" unless @ARGV;

my $coverage;
my @gff;
my $output = 'expressed_genes.tsv';
GetOptions(
	'c=s'=>\$coverage,
	'o=s'=>\$output,
	'g=s@{1,}'=>\@gff
);


## Creating coverage database
open COV, "<", "$coverage" or die "Can't open $coverage: $!\n";
my %coverage;
print "Working on coverage database...\n";
while (my $line = <COV>){
	chomp $line;
	if ($line =~ /^(\S+)\t(\d+)\t(\d+)/){
		my $contig = $1; #print "$contig\t";
		my $position = $2; #print "$position\t";
		my $cov = $3; #print "$cov\n";
		$coverage{$contig}{$position} = $cov;
	}
}

## Working through CDS annotations; disregarding tRNA and rRNAs...
print "Working on GFF files...\n";
open OUT, ">", "$output" or die "Can't create $output: $!\n";
print OUT "Protein\tProduct\tSum\tLength (in nt)\tAverage RNAseq depth\n";

my $sum;
my $size;
my $average;
my $cg;
my $locus_tag;
my $product;
my %products;
my %exons;
my $left;
my $right;
while (my $gff = shift@gff){
	open GFF, "<", "$gff" or die "Can't open $gff: $!\n";
	while (my $line = <GFF>){
		chomp $line;
		$sum = 0; $size = 0; $average = 0;
		if ($line =~ /^#/){ next; }
		elsif ($line =~ /^(\S+).*(Genbank|EMBL|RefSeq)\sgene.*locus_tag=(\w+)/){
			$cg = $1;
			$locus_tag = $3;
		}
		elsif ($line =~ /^(\S+).*(Genbank|EMBL|RefSeq)\sCDS\s(\d+)\s+(\d+).*;product=(.*?);/){
			$cg = $1;
			$left = $3;
			$right= $4;
			$product = $5;
			$products{$locus_tag}[0] = $product;
			$products{$locus_tag}[1] = $cg;
			push (@{$exons{$locus_tag}}, "$left\t$right");
		}
	}
}

for my $tag (keys %products){
	while (my $nt = shift@{$exons{$tag}}){
		my ($start, $end) = split ("\t", $nt);
		for my $pos ($start..$end){
			$sum += "$coverage{$products{$tag}[1]}{$pos}";
			$size++;
		}
	}
	$average = sprintf ("%.2f", $sum/$size); # unless $size == 0;
	print OUT "$tag\t$products{$tag}[0]\t$sum\t$size\t$average\n";
	$sum = 0; $size = 0; $average = 0;
}
