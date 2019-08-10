#!/usr/bin/perl

use strict; use warnings;

my $evalue = 1e-10;

open IN, "<Table_S2_v3_extended.txt";
open OUT, ">Table_S2_v3_extended.$evalue.txt";

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^A3770/){
		my @hits = split("\t", $line);
		for my $i (6,8,11,13,16,18,21,23,26,28,31,33,36,38,41,43,46,48,51,53,56,58,61,63,66,68){
			if ($hits[$i] >= $evalue) {$hits[$i]  = '---'; $hits[$i-1] = '---';}
		}
		for (my $x = 0; $x <= 67; $x++ ){print OUT "$hits[$x]\t";}
		print OUT "$hits[68]\n";
	}
	else {print OUT "$line\n";}
}