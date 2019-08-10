#!/usr/bin/perl

use strict;
use warnings;

my $usage = "USAGE = get_matrix.pl COCCO.vs.CCMP.txt COCCO.vs.CCMP.txt.metrics CCMP.list COCCOvsCCMP.matrix";
die "\n$usage\n\n" unless @ARGV;

open IN, "<$ARGV[0]";
open METRICS, "<$ARGV[1]";
open LIST, "<$ARGV[2]";
open OUT, ">$ARGV[3]";

## Creating matrix header
my %seen; my @reference;
while (my $line = <LIST>){
	chomp $line;
	if ($line =~ /^\S+\s+(\S+)/){
		if (!exists $seen{$1}){$seen{$1} = $1; push (@reference, $1);}
	}
}
my @refs = @reference;
print OUT "Contig";
while (my $ctg = shift@reference){print OUT ",$ctg";}
print OUT "\n";

## Creating matrix
my %hash; my $contig; my @contigs; my %hits;
while (my $line = <METRICS>){
	chomp $line;
	if ($line =~ /(\w+.?\w+).list.*blast hits = (\d+)/){
		$hits{$1} = $2;
	}
}	

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /(\w+.?\w+).list$/){
		$contig = $1;
		push (@contigs, $contig);
		$hash{$contig} = undef;
	}
	elsif ($line =~ /^(\S+)\s+(\d+)/){
		my $chromo = $1;
		my $number = $2;
		$hash{$contig}{$chromo} = $number;
	}
}
while (my $cg = shift@contigs){
	print OUT "$cg";
	foreach (@refs){
		if (exists $hash{$cg}{$_}){
			my $percent = ($hash{$cg}{$_}/$hits{$cg}*100);
			$percent = sprintf ("%.2f", $percent);
			print OUT ",$percent";}
		else {print OUT ",0";}
	}
	print OUT "\n";
}