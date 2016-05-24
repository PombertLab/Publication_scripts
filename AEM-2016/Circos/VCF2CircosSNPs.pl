#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'perl VCF2CircosSNPs.pl karyotype_accessions *.vcf';

my $db = shift@ARGV;
my %hash = ();
open HASH, "<$db";
while (my $tmp = <HASH>){
	chomp $tmp;
	if ($tmp =~ /^(\S+)\t(\S+)/){$hash{$1}=$2;}
}

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.circosSNPs";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^gi\|\d+\|\w+\|(\w+)\S+\|\s(\d+)/){
			my $acc = $1;
			my $pos = $2;
			print OUT "$hash{$acc} $pos $pos\n";
		}
	}
}