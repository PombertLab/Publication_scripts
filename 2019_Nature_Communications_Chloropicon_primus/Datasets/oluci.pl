#!/usr/bin/perl

open FA, "<Olucimarinus_231_v2.0.protein.fa" or die "Can't open Olucimarinus_231_v2.0.protein.fa";
open DEF, "<Olucimarinus_231_v2.0.defline.txt" or die "Can't open Olucimarinus_231_v2.0.defline.txt";
open OUT1, ">Olucimarinus.fasta";
open OUT2, ">Olucimarinus.products";

my %products;
while (my $line = <DEF>){
	if ($line =~ /^(\d+)\s+pdef\s+(.*)$/){$products{$1} = $2;}
}

while (my $line = <FA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){
		if (exists $products{$1}){
			print OUT1 ">$1\t$products{$1}\n";
			print OUT2 "$1\t$products{$1}\n";
		}
		else {
			print OUT1 ">$1\thypothetical protein\n";
			print OUT2 "$1\thypothetical protein\n";
		}
	}
	else {print OUT1 "$line\n";}
}


