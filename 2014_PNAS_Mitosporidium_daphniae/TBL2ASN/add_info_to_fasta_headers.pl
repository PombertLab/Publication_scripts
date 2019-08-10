#!/usr/bin/perl

use strict;
use warnings;

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.headers";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>.*$/){
			print OUT "$line [organism=Mitosporidium daphniae][strain=JFP-2104][lineage=cellular organisms; Eukaryota; Fungi/Metazoa group; Fungi; Microsporidia;][gcode=1]\n";
		}
		else {print OUT "$line\n";}
	}
}