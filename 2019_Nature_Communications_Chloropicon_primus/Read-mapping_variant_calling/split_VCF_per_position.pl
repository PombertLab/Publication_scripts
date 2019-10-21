#!/usr/bin/perl
## Version 0.1

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## SYNOPSIS: Splits VCF file per position on the contig. Useful to search for distored frequencies within a contig/chromosome.

die "\nUSAGE = sort_VCF_per_position.pl -start 1 -end 540000 -vcf *.vcf\n\n" unless @ARGV;

my $start;
my $end;
my @VCF;
GetOptions(
	'start=i' => \$start,
	'end=i' => \$end,
	'v|vcf=s@{1,}' => \@VCF
);

while (my $file = shift@VCF){
	open IN, "<$file"; $file =~ s/.vcf$//;
	open SUB, ">$file.${start}-${end}.vcf";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @columns = split("\t", $line);
			if (($columns[1]>=$start)&&($columns[1]<=$end)){print SUB "$line\n";}
		}
	}
}