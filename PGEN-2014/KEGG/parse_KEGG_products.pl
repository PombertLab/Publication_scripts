#!/usr/bin/perl

use strict;
use warnings;

my $usage = 'CORRECT USAGE = perl parse_product.pl KEGG_list File.parsed';

### Here we fill the product database
my $hash = shift@ARGV or die $usage;
my %products = ();
open HASH, "<$hash";
while (my $line = <HASH>){
	chomp $line;
	if ($line =~ /^(\d+)(\t.*)$/){
		my $key = $1;
		my $value = $2;
		$products{$key} = $value;
	}
}
### Here we work on the file
while (my $file = shift@ARGV){
	open IN, "<$file";
	open OUT, ">$file.KEGG_products";
	while (my $text = <IN>) {
		chomp $text;
		if ($text =~ /^(a\d+\;\d+)\t(\d+)\t(\S+)/){
			my $rna = $1;
			my $proteinid = $2;
			my $evalue = $3;
		
			if (exists $products{$proteinid}){
				print OUT "$rna\t$evalue\t$proteinid\t$products{$proteinid}\n";
			}
			else{
				print OUT "$rna\t$evalue\t$proteinid\tNo similar product found\n";
			}
		}
	}
}
