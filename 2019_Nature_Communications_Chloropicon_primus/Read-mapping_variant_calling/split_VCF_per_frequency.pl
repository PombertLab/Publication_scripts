#!/usr/bin/perl
# Pombert lab 2019
my $version = '0.1';
my $name = 'split_VCF_per_frequency.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Splits VCF file per allelic frequency. Useful to search for distored frequencies within a contig/chromosome.
USAGE		sort_VCF_per_frequency.pl -min 43 -max 57 -vcf *.vcf
OPTIONS
die "$options\n" unless @ARGV;

my $min = 43;
my $max = 57;
my @VCF;
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'v|vcf=s@{1,}' => \@VCF
);

while (my $file = shift@VCF){
	open IN, "<$file"; $file =~ s/.vcf$//;
	open N2, ">$file.$min.$max.N2.vcf";
	open N3, ">$file.$min.$max.N3.vcf";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @columns = split("\t", $line);
			my @stats = split(":", $columns[9]);
			$stats[6] =~ s/%//; $stats[6] = sprintf("%.03d", $stats[6]);
			if ($columns[3] eq 'N'){next;} ## Skipping masked nucleotides
			elsif (($stats[6] >= $min) && ($stats[6] <= $max)){print N2 "$line\n";}
			elsif (($stats[6] < $min) && ($stats[6] > 20)){print N3 "$line\n";}
			elsif (($stats[6] > $max) && ($stats[6] <= 80)){print N3 "$line\n";} 
		}
	}
}
