#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'get_products.pl filtered_blastp_output ## in table format';
die $usage unless @ARGV;

open HASH, "<uniprot.hash"; ## large file requires memory...
my %hash= ();

system "echo 'Filling the UniProtDB, this may take a while...'";
while (my $dbkey = <HASH>){
	chomp $dbkey;
	if ($dbkey =~ /(\S+)\t(.*$)/){
		my $key = $1;
		my $desc = $2;
		$hash{$key}=$desc;
	}
}

while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.products";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(DI09_\d+p\d+)\ttr\|(\S+)\|/){
			my $protein = $1;
			my $uniprot = $2;
			print OUT "$protein\t$hash{$uniprot}\n";
		}
	}
}