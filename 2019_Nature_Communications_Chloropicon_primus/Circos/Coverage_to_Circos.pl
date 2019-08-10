#!/usr/bin/perl
## Pombert Lab, 2019
my $name = 'Coverage_to_Circos.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Generates a sequencing depth distribution using sliding windows for Circos
EXAMPLE		Coverage_to_Circos.pl -f *.coverage -s 500 -w 1000

-h (-help)	Displays list of options	
-f (--file)	Coverage input file(s) ## Samtools depth -aa function
-o (--ouput)	Output file names prefix [Default: genome]
-c (--color)	Color for genotype [Default: black]
-s (--step)	Size of the steps between windows [Default: 500]
-w (--window)	Width of the sliding windows [Default: 1000]
OPTIONS
die "$options\n" unless @ARGV;

my $help;
my @file;
my $output = 'genome';
my $color = 'black';
my $slide = 500; my $window = 1000;

GetOptions(
	'h|help' => \$help,
	'f|file=s@{1,}' => \@file,
	'o|ouput=s' => \$output,
	'c|color=s' => \$color,
	's|slide=i' => \$slide,
	'w|window=i' => \$window
);

die $options if $help;

my %sequences;
my %seq_list;
my @list;
## Working on coverage file
while (my $file = shift @file) {
	open IN, "<$file" or die "Can\'t open $file";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\S+)\t(\d+)\t(\d+)/){
			my $chr = $1;
			my $pos = $2;
			my $cov = $3;
			$sequences{$chr}{$pos} = $cov;
			unless (exists $seq_list{$chr}){push (@list, $chr);}
			$seq_list{$chr} += 1;
		}
	}
	close IN;
}
## Creating a "karyotype" file for Circos
open COV, ">$output.cov"; print COV '#chr START END Coverage'."\n";
open KAR, ">$output.genotype"; print KAR '#chr - ID LABEL START END COLOR'."\n";
my $id = 0; my $mmax; my $mmin;

## Working on sequences
while (my $sequence = shift @list){
	my $len = $seq_list{$sequence};
	print "\nWorking on $sequence ... \n";
	print "$sequence size = $len\n"; ## Debugging
	my $terminus = $len - 1;
	my $max; my $min;
	$id ++;
	print KAR "chr - $sequence $id 0 $terminus $color\n";
	for(my $i = 1; $i < $len - $window; $i += $slide){
		my $sum_cov = 0; my $pos = 0;
		for ($i..$i+$window-1){
				$sum_cov += $sequences{$sequence}{$_};
				$pos++;
		}
		my $av_cov = sprintf("%.2f", $sum_cov/$pos);
		my $end = $i + $pos - 1;
		print COV "$sequence $i $end $av_cov\n";
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
print "\nMin Cov = $mmin \n";
print "Max Cov = $mmax \n";
close COV;
close KAR;
exit;
