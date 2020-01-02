#!/usr/bin/perl
## Pombert Lab, 2017
my $name = 'seq_depth.pl';
my $version = 0.1;

use strict; use warnings;

## Defining options
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Compares the sequencing depths of contigs to that their assembly average
		The sequencing depth information is obtained from samtools: samtools depth -aa file.bam > output.coverage
USAGE		seq_depth.pl *.coverage
NOTE		seq_depth.pl is now integrated in get_SNPs.pl (https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl)

OPTIONS
die "$usage" unless @ARGV;

## Working on files
while (my $file = shift@ARGV){
	open COV, "<$file";
	$file =~ s/.coverage$//;
	open DEPTH, ">$file.depth";
	my $total = 0; my $covered = 0; my $nocov = 0; my $max = 0; my $sumcov; my $sn =0;
	print DEPTH "Contig\tAverage depth\tAverage (all contigs)\tDifference\tRatio Contig\/Average\n";
	my %data = (); my @contigs = ();
	while (my $line = <COV>){
		chomp $line;
		$total++; 
		if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)/){
			my $contig = $1; my $position = $2; my $coverage = $3;
			$sumcov += $coverage;
			if ($coverage >= 1) {
				$covered++;
				if ($coverage > $max){$max = $coverage;}
			}
			else {$nocov++;}
			if (exists $data{$contig}){
				$data{$contig}[0] += 1; ## Keeping track of contig size 
				$data{$contig}[1] += $coverage; ## Keeping track of cumulative sequencing depth 
			}
			else{
				$data{$contig}[0] = 1; ## Initializing new contig 
				$data{$contig}[1] += $coverage; ## Keeping track of cumulative sequencing depth 
				push (@contigs, $contig);
			}
		}
	}
	my $avc = sprintf("%.2f", ($sumcov/$total));
	while (my $tmp = shift@contigs){ ## Printing sequencing depths per contig
		my $average = ($data{$tmp}[1]/$data{$tmp}[0]);
		$average = sprintf("%.2f", $average);
		my $diff = $average - $avc; $diff = sprintf("%.2f", $diff);
		my $ratio = $average/$avc; $ratio = sprintf("%.2f", $ratio);
		print DEPTH "$tmp\t$average\t$avc\t$diff\t$ratio\n";
	} 
}