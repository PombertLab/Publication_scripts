#!/usr/bin/perl

open IN, "<$ARGV[0]";
open OUT, ">Otauri.products";
while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^>(\S+)\s(.*)\s\[Ostreococcus tauri\]/){
		print OUT "$1\t$2\n";
	}
}