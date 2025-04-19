#!/usr/bin/env perl
## Pombert Lab, 2019

my $name = 'Coverage_to_Circos.pl';
my $version = 0.2;
my $updated = '2023-10-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Generates a sequencing depth distribution using sliding windows for Circos

EXAMPLE		${name} \\
		  -f *.coverage \\
		  -n \\
		  -s 500 \\
		  -w 1000

-f (--file)	Coverage input file(s) ## Samtools depth -aa function
-o (--output)	Output file names prefix [Default: genome]
-n (--normal)	Normalize by the average sequencing depth to calculate # of alleles
-p (--ploidy)	Ploidy to use for normalization [Default: 1]
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]
-g (--gen)	Produce a genotype/karyotype file for Circos
-c (--color)	Color for genotype [Default: black]
OPTIONS
die "\n$options\n" unless @ARGV;

my @file;
my $output = 'genome';
my $normalized;
my $ploidy = 1;
my $slide = 500;
my $window = 1000;
my $genotype;
my $color = 'black';

GetOptions(
	'f|file=s@{1,}' => \@file,
	'o|output=s' => \$output,
	'n|normal' => \$normalized,
	'p|ploidy=i' => \$ploidy,
	's|slide=i' => \$slide,
	'w|window=i' => \$window,
	'g|gen' => \$genotype,
	'c|color=s' => \$color
);

my %sequences;
my %seq_lengths;
my %metrics;

## Working on coverage file
while (my $file = shift @file) {

	open IN, "<", $file or die "Can\'t open $file";

	while (my $line = <IN>){

		chomp $line;

		if ($line =~ /^#/){
			next;
		}

		else {

			my @data = split("\t", $line);
			my $chr = $data[0];
			my $pos = $data[1];
			my $cov = $data[2];

			$sequences{$file}{$chr}{$pos} = $cov;
			$seq_lengths{$file}{$chr}{'len'}++;

			$metrics{$file}{'sumlen'}++;
			$metrics{$file}{'sumcov'} += $cov;

		}
	}
	close IN;

}

## Creating a coverage file for Circos

my $cov_file = $output.'.cov';
open COV, '>', $cov_file or die "Can't create $cov_file: $!\n";
print COV '#chr START END Coverage'."\n";

if ($genotype){
	my $gen_file = $output.'genotype';
	open KAR, '>', $gen_file or die "Can't create $gen_file: $!\n";
	print KAR '#chr - ID LABEL START END COLOR'."\n";
}

my $id = 0;
my $mmax;
my $mmin;

## Working on sequences

foreach my $file (sort (keys %sequences)){

	foreach my $sequence (sort(keys %{$sequences{$file}})){

		my $len = $seq_lengths{$file}{$sequence}{'len'};

		my $sumlen = $metrics{$file}{'sumlen'};
		my $sumcov = $metrics{$file}{'sumcov'};

		my $average_depth = $sumcov/$sumlen;

		print "\nWorking on $sequence ... \n";
		print "$sequence size = $len\n"; ## Debugging
		
		my $terminus = $len - 1;
		my $max;
		my $min;
		$id ++;

		if ($genotype){
			print KAR "chr - $sequence $id 0 $terminus $color\n";
		}

		for (my $i = 1; $i < $len - $window; $i += $slide){

			my $sum_cov = 0;
			my $pos = 0;

			for ($i..$i+$window-1){
					$sum_cov += $sequences{$file}{$sequence}{$_};
					$pos++;
			}

			my $av_cov = sprintf("%.2f", $sum_cov/$pos);
			my $normalized_cov = sprintf("%.2f", ($av_cov/$average_depth) * $ploidy);

			my $end = $i + $pos - 1;

			if ($normalized){
				print COV "$sequence $i $end $normalized_cov\n";
			}
			else{
				print COV "$sequence $i $end $av_cov\n";
			}

			if (!defined $max){$max = $av_cov;}
			if (!defined $mmax){$mmax = $av_cov;}
			if (!defined $min){$min = $av_cov;}
			if (!defined $mmin){$mmin = $av_cov;}

			if ($av_cov > $max){$max = $av_cov;}
			if ($av_cov > $mmax){$mmax = $av_cov;}
			if ($av_cov < $min){$min = $av_cov;}
			if ($av_cov < $mmin){$mmin = $av_cov;}

		}

		print "Min/Max Cov on $sequence = $min / $max \n";

	}

}

print "\nMin Cov = $mmin \n";
print "Max Cov = $mmax \n";

close COV;
if ($genotype){
	close KAR;
}

exit;