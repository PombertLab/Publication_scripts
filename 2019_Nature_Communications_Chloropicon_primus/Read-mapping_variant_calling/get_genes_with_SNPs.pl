#!/usr/bin/perl
## Pombert Lab, 2019
my $version = 0.2;
my $name = 'get_genes_with_SNPs.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a tab-delimited list of variants and their location (CDS, RNA, intron, intergenic).
USAGE		get_genes_with_SNPs.pl -vcf *.vcf -embl *.embl  -tbl *.tbl -p Verified_products_ALL.list -o SNP.tsv 

OPTIONS:
-vcf	Variants identified in VCF format
-embl	CDS annotations in EMBL format ## Converted from WebAppolo GFF files, which do not contain RNAs
-tbl	rRNA/tRNA annotations in TBL format
-p	List of locus tags and their products
-o	Desired output file name
OPTIONS
die "$usage\n" unless @ARGV; 

my @embl;
my @vcf;
my @tbl;
my $products = 'Verified_products_ALL.list';
my $out = 'SNP.tsv';
GetOptions(
	'embl=s@{1,}' => \@embl,
	'vcf=s@{1,}' => \@vcf,
	'tbl=s@{1,}' => \@tbl,
	'p=s' => \$products,
	'o=s' => \$out
);

## Creating database of locus_tags => products
my %prod;
open IN, "<$products";
open OUT, ">$out";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^(\S+)\t(.*)$/){$prod{$1} = $2;}
}

## Creating databases of CDS and introns
my %hashes; my %introns; my $key; my $tag;
while (my $file = shift @embl){
	open EMBL, "<$file";
	$file =~ s/.embl$//;
	while (my $line = <EMBL>){
		chomp $line;
		if ($line =~ /^FT\s+gene/){$key = 'gene';} ## EMBL
		elsif ($line =~ /^FT\s+\/locus_tag=\"(.*)\"/){$tag = $1;} ## EMBL
		elsif ($line =~ /^FT\s+intron\s+(complement\()?(\d+)\.\.(\d+)/){
			my $start = $2; my $end = $3;
			for ($start..$end){
				push(@{$introns{$file}{$_}},  $tag);
			}
		}
		elsif ($line =~ /^FT\s+CDS\s+(complement\()?(join\()?([0-9.,]+)/){
			my @posits = split(",", $3);
			while (my $tmp = shift@posits){
				if ($tmp =~ /(\d+)\.\.(\d+)/){
					my $start = $1; my $end = $2;
					for ($start..$end){
						push(@{$hashes{$file}{$_}},  $tag);
					}
				}
			}
		}
	}
	my @CDS = keys %{ $hashes{$file} }; my $cds_count = scalar @CDS;
	my @INTRONS =  keys %{ $introns{$file} }; my $introns_count = scalar @INTRONS;
	print "CDS (exon-only) bases in $file\t$cds_count\n";
	print "Intron bases in $file\t$introns_count\n";
}
## Creating tRNA and rRNA databases
my %rRNAs; my %tRNAs;
while (my $file = shift @tbl){
	open TBL, "<$file";
	$file =~ s/.tbl$//;
	while (my $line = <TBL>){
		chomp $line;
		if ($line =~ /^\s+locus_tag\s+(A3770_\w+)/){$tag = $1;} ## EMBL
		elsif ($line =~ /^(\d+)\s+(\d+)\s+rRNA/){
			my $start = $1; my $end = $2;
			for ($start..$end){
				push(@{$rRNAs{$file}{$_}},  $tag);
			}
		}
		elsif ($line =~ /^(\d+)\s+(\d+)\s+tRNA/){
			my $start = $1; my $end = $2;
			for ($start..$end){
				push(@{$tRNAs{$file}{$_}},  $tag);
			}
		}
	}
	my @rRNA = keys %{ $rRNAs{$file} }; my $rRNAs_count = scalar @rRNA;
	my @tRNA =  keys %{ $tRNAs{$file} }; my $tRNAs_count = scalar @tRNA;
	print "rRNA bases in $file\t$tRNAs_count\n";
	print "tRNA bases in $file\t$rRNAs_count\n";
}
## Parsing VCF file(s)
print OUT "Contig\tPosition\tReference\tAlternative base\tAllelic Frequency\tType\tLocus_tag\tProduct\n";
while (my $vcf = shift@vcf){
	open VCF, "<$vcf";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @vcf = split("\t", $line); my @params = split(":", $vcf[9]); my $freq = $params[6]; my $evalue = $params[7];
			my $contig = $vcf[0]; my $position = $vcf[1]; my $ref = $vcf[3]; my $alt = $vcf[4];
			if (exists $hashes{$contig}{$position}){
				print OUT "$contig\t$position\t$ref\t$alt\t$freq\tCDS\t$hashes{$contig}{$position}[0]\t$prod{$hashes{$contig}{$position}[0]}\n";
			}
			elsif (exists $introns{$contig}{$position}){
				print OUT "$contig\t$position\t$ref\t$alt\t$freq\tintron\t$introns{$contig}{$position}[0]\t$prod{$introns{$contig}{$position}[0]}\n";
			}
			elsif (exists $tRNAs{$contig}{$position}){
				print OUT "$contig\t$position\t$ref\t$alt\t$freq\tintron\t$tRNAs{$contig}{$position}[0]\t$prod{$tRNAs{$contig}{$position}[0]}\n";
			}
			elsif (exists $rRNAs{$contig}{$position}){
				print OUT "$contig\t$position\t$ref\t$alt\t$freq\tintron\t$rRNAs{$contig}{$position}[0]\t$prod{$rRNAs{$contig}{$position}[0]}\n";
			}
			else{
				print OUT "$contig\t$position \t$ref\t$alt\t$freq\tintergenic\tN\/A\tN\/A\n";
			}
		}
	}
}
