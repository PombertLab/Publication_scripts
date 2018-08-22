#!/usr/bin/perl

use strict;
use warnings;

my $usage = "
USAGE = script.pl evalue
\n";
die "$usage" unless@ARGV;
my $evalue = shift@ARGV;
my $scifmt = sprintf("%.100f", $evalue);

open IN, "<list.txt";
open IN2, "<prasinos.hits";
open IN3, "<others.hits";
open IN4, "<Verified_products_ALL.list";
open OUT, ">In_Prasinos.$evalue.txt";
open OUT2, ">Not_in_Prasinos_but_in_other_algae.$evalue.txt";
open OUT3, ">Uniques_to_CCMP.$evalue.txt";

my %prasinos = ();
my %others = ();
my %products = ();

while (my $line = <IN4>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$products{$1}=$2;}
}

while (my $line = <IN2>){
	chomp $line;
	if ($line =~ /^(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/){
		my $locus = $1; my $hit = $2; my $eval = $3;
		my $ev = sprintf("%.100f", $eval);
		if ($ev <= $scifmt){
			$prasinos{$locus} = $locus;
		}
	}
}
while (my $line = <IN3>){
	chomp $line;
	if ($line =~ /^(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)/){
		my $locus = $1; my $hit = $2; my $eval = $3;
		my $ev = sprintf("%.100f", $eval);
		if ($ev <= $scifmt){
			$others{$locus} = $locus;
		}
	}
}
while (my $line = <IN>){
	chomp $line;
	if (exists $prasinos{$line}){
		print OUT "$line\t$products{$line}\n";
	}
	elsif (exists $others{$line}){
		print OUT2 "$line\t$products{$line}\n";
	}
	else{
		print OUT3 "$line\t$products{$line}\n";
	}
}