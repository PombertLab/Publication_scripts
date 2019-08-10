#!/usr/bin/perl

use strict;
use warnings;

my $usage = "perl parse_BLASTP.pl min_evalue *.blastp";
die "$usage\n" unless @ARGV;

my $evalue = shift@ARGV;
my $scifmt = sprintf("%.100f", $evalue);

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.parsed";
	
	my %hits = ();
	my %evalues = ();
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/){
			my $locus = $1; my $hit = $2; my $eval = $3;
			if (exists $hits{$locus}){next;}
			else{
				$hits{$locus} = $hit;
				$evalues{$locus} = $eval;
			}
		}
	}
	open IN2, "<list.txt";
	while (my $line = <IN2>){
		chomp $line;
		if (exists $hits{$line}){
			my $ev = sprintf("%.100f", $evalues{$line});
			if ($ev <= $scifmt){
				print OUT "$line\t$hits{$line}\t$evalues{$line}\n";
			}
			else{
				print OUT "$line\t---\t---\n";
			}
		}
		else{
			print OUT "$line\t---\t---\n";
		}
	}
}
	