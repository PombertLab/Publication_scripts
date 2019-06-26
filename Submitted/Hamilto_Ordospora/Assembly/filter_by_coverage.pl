#!/usr/bin/perl
## Pombert Lab, 2017
my $name = 'filter_by_coverage.pl';
my $version = 0.1;

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Parses contigs by their average sequencing depth
USAGE		filter_by_coverage.pl -min 20 -max 200 -fa *.fasta -cov *.coverage -out depth_summary.txt

OPTIONS:
-min		Desired minimum sequencing depth [default: 20]
-max		Desired maximum sequencing depth [default: 200]
-fa		Fasta file(s) to parse
-cov		Coverage files obtained by read mapping with get_SNPs.pl
-out		Depth summary output [default: depth_summary.txt]
OPTIONS
die "$options\n" unless @ARGV;

my $min = 20;
my $max = 200;
my @fasta;
my @cov;
my $output = 'depth_summary.txt';
GetOptions(
	'min=i' => \$min,
	'max=i' => \$max,
	'fa|fasta=s@{1,}' => \@fasta,
	'cov|coverage=s@{1,}' => \@cov,
	'out=s' => \$output
);

my %coverage;
while (my $coverage = shift@cov){
	open IN, "<$coverage";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\S+)\t(\d+)\t(\d+)/){
			my $locus = $1;
			my $position = $2;
			my $depth = $3;
			$coverage{$locus}[0] = $position; 
			$coverage{$locus}[1] += $depth;
		}
	} 
}

open OUT, ">$output";
while (my $fasta = shift@fasta){
	open IN, "<$fasta";
	$fasta =~ s/\.w+//;
	open MIN, ">$fasta.${min}-X.fasta";
	open MED, ">$fasta.${min}_to_${max}X.fasta";
	open MAX, ">$fasta.${max}+X.fasta";
	my %sequences;
	my $contig;
	my @contigs;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)(.*)$/){
			push (@contigs, $1);
			$contig = $1;
		}
		else{
			$sequences{$contig} .= $line;
		}
	}
	while (my $cg = shift@contigs){
		my $avcov = $coverage{$cg}[1]/$coverage{$cg}[0];
		if ($avcov < $min){
			print OUT "UNDER\t$cg\tlength = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MIN ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MIN "$sq\n";}
		}
		elsif ($avcov > $max){
			print OUT "OVER\t$cg\tlength = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MAX ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MAX "$sq\n";}
		}
		else{
			print OUT "BETWEEN\t$cg\tlength = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MED ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MED "$sq\n";}
		}
	}
}
close OUT;
open IN, "<$output";
my $under = 0; my $between = 0; my $over = 0;
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^UNDER.*length = (\d+)/){$under+=$1;}
	elsif ($line =~ /^BETWEEN.*length = (\d+)/){$between+=$1;}
	elsif ($line =~ /^OVER.*length = (\d+)/){$over+=$1;}
}

my $sum = $under + $between + $over;
print "\nUnder ${min}X\t\t$under bp\n";
print "In-between ${min}-${max}X\t$between bp\n";
print "Over ${max}X\t\t$over bp\n";
print "Total\t\t\t$sum bp\n\n";
