#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = filter_by_coverage.pl -min min -max max -fa *.fasta -cov *.coverage";
die "\n$usage\n\n" unless @ARGV;

my $min;
my $max;
my @fasta;
my @cov;

GetOptions(
	'min=s' => \$min,
	'max=s' => \$max,
	'fa|fasta=s@{1,}' => \@fasta,
	'cov|coverage=s@{1,}' => \@cov
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

while (my $fasta = shift@fasta){
	open IN, "<$fasta";
	open MIN, ">$fasta.mincov";
	open MED, ">$fasta.medcov";
	open MAX, ">$fasta.maxcov";
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
			print "UNDER\t$cg\tlenght = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MIN ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MIN "$sq\n";}
		}
		elsif ($avcov > $max){
			print "OVER\t$cg\tlenght = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MAX ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MAX "$sq\n";}
		}
		else{
			print "BETWEEN\t$cg\tlenght = $coverage{$cg}[0]\tAverage sequencing depth = $avcov\n";
			print MED ">$cg\n";
			my @seq = unpack("(A60)*", $sequences{$cg});
			while (my $sq = shift@seq){print MED "$sq\n";}
		}
	}
}
	
