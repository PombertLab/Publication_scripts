#!/usr/bin/perl

## Reformats output from KEGG Orthology database [http://www.genome.jp/kegg/ko.html, Get title option] into tab-delimited lists

use strict;
use warnings;

while (my $file = shift @ARGV){
	open IN, "<$file";
	open OUT, ">$file.out";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^ko\:(\w+).*;\s(.*)\s\[EC:(.*)\]/){
			my $ko = $1;
			my $product = $2;
			my $ec = $3;
			print OUT "$ko\t$ec\t$product\n";
		}
		elsif ($line =~ /^ko\:(\w+).*;\s(.*)/){
			my $ko = $1;
			my $product = $2;
			my $ec = undef;
			print OUT "$ko\t---\t$product\n";
		}
	}
}
close IN;
close OUT;
exit;