#!/usr/bin/perl

open IN, "<$ARGV[0]";
open OUT, ">hash.formatted";

my $memory = 'ABC';

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^(K\d{5})\t(.*)$/){
		my $ko = $1;
		my $protein = $2;
		if ($ko eq $memory){
			print OUT ", $protein";
		}
		elsif($ko ne $memory){
			print OUT "\n$ko\tCOCSUDRAFT_$protein";
			$memory = $ko;
		}
	 }
}
print OUT "\n";
	