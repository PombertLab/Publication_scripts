#!/usr/bin/perl

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = 'USAGE = get_genes_with_SNPs.pl -embl *.embl -vcf *.vcf > SNP.tsv';
die "\n$usage\n\n" unless @ARGV; 

my @embl; my @vcf; my $products = 'Verified_products_ALL.list';
GetOptions(
	'embl=s@{1,}' => \@embl,
	'vcf=s@{1,}' => \@vcf,
	'p=s' => \$products
);

my %prod;
open IN, "<$products";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^(\S+)\t(.*)$/){$prod{$1} = $2;}
}

my %hashes; my $key; my $tag;
while (my $file = shift @embl){
	open EMBL, "<$file";
	$file =~ s/.embl$//;
	while (my $line = <EMBL>){
		chomp $line;
		if ($line =~ /^FT\s+gene/){$key = 'gene';} ## EMBL
		elsif ($line =~ /^FT\s+\/locus_tag=\"(.*)\"/){$tag = $1;} ## EMBL
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
}
print "Contig\tPosition\tReference\tAlternative base\tAllelic Frequency\tLocus_tag\tProduct\n";
while (my $vcf = shift@vcf){
	open VCF, "<$vcf";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @vcf = split("\t", $line); my @params = split(":", $vcf[9]); my $freq = $params[6]; my $evalue = $params[7];
			my $contig = $vcf[0]; my $position = $vcf[1]; my $ref = $vcf[3]; my $alt = $vcf[4];
			if (exists $hashes{$contig}{$position}){
				print "$contig\t$position\t$ref\t$alt\t$freq\t$hashes{$contig}{$position}[0]\t$prod{$hashes{$contig}{$position}[0]}\n";
			}
			else{
				print "$contig\t$position \t$ref\t$alt\t$freq\tintergenic\tN\/A\n";
			}
		}
	}
}