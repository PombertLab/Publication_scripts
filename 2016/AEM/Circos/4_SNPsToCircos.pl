#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'SNPsToCircos.pl sliding_window_size circos_chr_name *.windows';
die $usage unless @ARGV;

my $win = shift@ARGV;
my $name = shift@ARGV;
my $chr = 0;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.SNPs";
	$chr++;
	my $start = 0;
	my $end = (-1);
	my $count = 0;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\d+)/){
			$count++;
			if ($count == $win){
				$end += $win;
				print OUT "$name$chr $start $end $line\n";
				$start += $win;
				$count = 0;
			}
			else{$end += $win;
			print OUT "$name$chr $start $end $line\n";
		}
		}
	}
}