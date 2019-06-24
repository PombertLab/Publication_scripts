#!/usr/bin/perl

use strict;
use warnings;

die "\nUSAGE = seq_depth.pl *.coverage\n" unless @ARGV;

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