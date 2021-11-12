#!/usr/bin/perl
## Pombert Lab, Illinois Tech, 2021

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $script_name = 'gesamt_table.pl';
my $version = '0.1';

my $USAGE=<<"OPTIONS";

NAME		${script_name}
VERSION		${version}
SYNOPSIS	Generate locus_tag and QScore table from gesamt results with predicted function of protein from ncbi cds fasta

EXAMPLE		${script_name} \\
		  -i *.gesamt \\
		  -n cds.fasta \\
		  -o table.tsv \\
		  -qs 0.3

OPTIONS:
-i (--in)	Gesamt files in the order of columns in table
-n (--ncbi)	CDS file from NCBI in fasta format
-o (--out)	Output name (table in .tsf format)
-qs (--QScore) 	Minimum quality score to filter gesamt files

OPTIONS
die "$USAGE" unless @ARGV;

my @files;
my $ncbi;
my $out;
my $minQS = 0.3;

GetOptions(
	'in|i=s{1,}' => \@files,
	'n|ncbi=s' => \$ncbi,
	'out|o=s' => \$out,
	'QScore|qs=f{1}' => \$minQS,
);

## Grabbing products from NCBI cds file
open NCBI, "<", "$ncbi" or die "Can't open $ncbi: $!\n";

my %products;
my $sequence;
my @sequence_loci;
my $sequence_length;

while (my $line = <NCBI>){

	chomp $line;

	## FASTA header with defined protein product
	if ($line =~ /^>.*gene=(\w+).*protein=(.*?)\]/){
		my $sequence_locus = $1;
		my $product = $2;
		$sequence_locus =~ s/i$//;
		$products{$sequence_locus}{'product'} = $product;
		push (@sequence_loci, $sequence_locus);

		if ($sequence){
			$sequence_length = (length $sequence)/3 - 1; ## -1 to remove stop codon
			$products{$sequence_loci[-2]}{'length'} = $sequence_length;
		}
		$sequence = undef;
	}

	## FASTA header without defined protein product
	elsif ($line =~ /^>.*gene=(\w+)/){
		my $sequence_locus = $1;
		my $product = 'undefined product';
		$sequence_locus =~ s/i$//;
		$products{$sequence_locus}{'product'} = $product;
		push (@sequence_loci, $sequence_locus);

		if ($sequence){
			$sequence_length = (length $sequence)/3 - 1; ## -1 to remove stop codon
			$products{$sequence_loci[-2]}{'length'} = $sequence_length;
		}
		$sequence = undef;
	}

	## Grabbing sequence
	else { $sequence .= $line; }
}
## Working on last sequence
$sequence_length = (length $sequence)/3 - 1;
$products{$sequence_loci[-1]}{'length'} = $sequence_length;

## Parsing GESAMT files
my %loci; 
open OUT, ">", "$out.tsv";

foreach my $file (@files){

	open IN, "<", "$file";

	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){ next; }
		my @columns = split(/\s+/, $line);
		my $QScore = $columns[3];
		my ($locus_tag) = $columns[8] =~ /(ECU\d+_\d+)/;

		## Parsing by QScore
		if ($QScore >= $minQS) {
			$loci{$locus_tag}{$file} = $QScore;
		}
	}
}

## Creating headers
print OUT "Locus_tag\tpredicted function\tlength (aa)";
foreach my $file (@files){ 
	my $header = $file;
	$header =~ s/\.(\w+)\.gesamt$//;
	print OUT "\t$header"; 
}
print OUT "\n";

## Working on table data content
foreach my $locus (sort (keys %loci)){
	print OUT "$locus";

	if (exists $products{$locus}{'product'}){ 
		print OUT "\t$products{$locus}{'product'}"; 
		print OUT "\t$products{$locus}{'length'}"; 
	}
	else { 
		print OUT "\thypothetical protein";
	}

	foreach my $file (@files){
		if (exists $loci{$locus}{$file}){
			print OUT "\t$loci{$locus}{$file}";
		}
		else { print OUT "\t---"; }
	}
	print OUT "\n";
}
