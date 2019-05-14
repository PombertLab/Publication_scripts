#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.2';
my $name = 'sort_SNPs.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Parse VarScan2 VCF allelic distributions for easy plotting with R or Excel
USAGE		sort_SNPs.pl -min 10 -max 90 -vcf *.vcf

OPTIONS:
-min	Mininum variant allelic frequency to keep [Default: 10]
-max	Maximum variant allelic frequency to keep [Default: 90]
-vcf	VarScan2 VCF files to parse
OPTIONS
die "$options\n" unless @ARGV;

my $min = 10;
my $max = 90;
my @VCF;
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'v|vcf=s@{1,}' => \@VCF
);

while (my $file = shift@VCF){
	open IN, "<$file"; $file =~ s/.vcf$//;
	open OUT1, ">$file.sorted.$min.$max.vcf";
	open OUT2, ">$file.sorted.$min.$max.tsv";
	open OUT3, ">$file.unsorted.$min.$max.tsv"; 	## File to graph with excel or R
	open OUT4, ">$file.gaussian.$min.$max.tsv"; 	## Data for Gaussian plot 
	my @percents; my %percents;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @columns = split("\t", $line);
			my @stats = split(":", $columns[9]);
			$stats[6] =~ s/%//; $stats[6] = sprintf("%.03d", $stats[6]);
			if ($columns[3] eq 'N'){next;} ## Skipping masked nucleotides
			elsif (($stats[6] >= $min) && ($stats[6] <= $max)){
				my $alternate = 100 - $stats[6]; ## Calculating frequency for alternate nucleotide, works only for diploid SNPs
				push (@percents, $stats[6]);
				push (@percents, $alternate);
				$percents{$stats[6]} += 1;
				$percents{$alternate} += 1;
				print OUT1 "$line\n";
			} 
		}
	}
	my @unsort = @percents;
	@percents = sort@percents;
	@percents = reverse@percents;
	while (my $pc = shift@percents){print OUT2 "$pc\n";}
	while (my $un = shift@unsort){print OUT3 "$un\n";}
	for (keys %percents){print OUT4 "$_\t$percents{$_}\n";}
}